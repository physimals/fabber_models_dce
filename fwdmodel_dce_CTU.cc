/*  fwdmodel_dce_ctu_NLLS.cc - Implementation of the compartment tissue update model

https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.26324

 Moss Zhao - IBME, Oxford

 Copyright (C) 2018 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce_CTU.h"

#include "fabber_core/easylog.h"
#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include <math.h>

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DCE_CTU_FwdModel> DCE_CTU_FwdModel::registration("dce_CTU");

static OptionSpec OPTIONS[] = {
    { "delt", OPT_FLOAT, "Time resolution between volumes, in minutes", OPT_REQ, "" },
    { "fa", OPT_FLOAT, "Flip angle in degrees.", OPT_REQ, "" },
    { "tr", OPT_FLOAT, "Repetition time (TR) In seconds.", OPT_REQ, "" },
    { "r1", OPT_FLOAT, "Relaxivity of contrast agent, In s^-1 mM^-1.", OPT_REQ, "" },

    { "Fp", OPT_FLOAT, "Flow in min-1", OPT_NONREQ, "0.3" },
    { "Vp", OPT_FLOAT, "Plasma volume in decimal between zero and one", OPT_NONREQ, "0.3" },
    { "Ve", OPT_FLOAT, "Extracellular space volume in decimal between zero and one", OPT_NONREQ, "0.3" },

    { "aif_type", OPT_STR, "Type of AIF must be aif_type=signal: User-supplied vascular MRI signal or aif_type=conc:User-supplied concentration curve", OPT_REQ, "none" },
    { "sig0", OPT_FLOAT, "Baseline signal. May be inferred.", OPT_NONREQ, "1000" },
    { "t10", OPT_FLOAT, "Baseline T1 value in seconds. May be inferred.", OPT_NONREQ, "1" },
    { "delay", OPT_FLOAT, "Delay time offset relative to AIF in minutes. May be inferred.", OPT_NONREQ, "0" },
    { "infer-delay", OPT_BOOL, "Infer the delay parameter", OPT_NONREQ, "" },
    { "infer-sig0", OPT_BOOL, "Infer baseline signal", OPT_NONREQ, "" },
    { "infer-T10", OPT_BOOL, "Infer T10 value", OPT_NONREQ, "" },
    { "convolution-method", OPT_STR, "Method to compute convolution, normal or iterative. Default is iterative", OPT_REQ, "" },
    { "use-log-params", OPT_BOOL, "Infer log values of parameters", OPT_NONREQ, "" },
    { "aif-file", OPT_FILE,
        "File containing single-column ASCII data defining the AIF. For aif=signal, this is the vascular signal curve. For aif=conc, it should be the blood plasma concentration curve",
        OPT_NONREQ, "none" },
    { "aif-hct", OPT_FLOAT, "Haematocrit value to use when converting an AIF signal to concentration. Used when aif=signal", OPT_NONREQ, "0.45" },
    { "aif-t1b", OPT_FLOAT, "Blood T1 value to use when converting an AIF signal to concentration. Used when aif=signal", OPT_NONREQ, "1.4" },
    { "" },
};

void DCE_CTU_FwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string DCE_CTU_FwdModel::ModelVersion() const
{
    string version = "Fabber DCE models: ";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void DCE_CTU_FwdModel::Initialize(FabberRunData &args)
{
    //cerr << "Initialize" << endl;
    // Mandatory parameters
    m_dt = args.GetDouble("delt", 0);
    m_FA = args.GetDouble("fa", 0, 90) * M_PI / 180; // convert from degree to radians
    m_TR = args.GetDouble("tr");
    m_r1 = args.GetDouble("r1");
    m_aif_type = args.GetString("aif_type");

    // Optional parameters
    initial_Fp = args.GetDoubleDefault("Fp", 0.3);
    initial_PS = args.GetDoubleDefault("PS", 0.3);
    initial_Vp = args.GetDoubleDefault("Vp", 0.3);

    // Case where we need to convert MRI signal into Gd concentration
    // Appendix A in https://www.sciencedirect.com/science/article/pii/S0730725X10001748
    if (m_aif_type == "signal")
    {
        ColumnVector aif_sig = read_ascii_matrix(args.GetString("aif-file"));
        // Convert AIF to concentration curve
        double aif_t1 = args.GetDoubleDefault("aif-t1", 1.4);
        double aif_hct = args.GetDoubleDefault("aif-hct", 0.45);
        double aif_sig0 = aif_sig(1);
        m_aif = aif_sig;
        for (int t = 0; t < m_aif.Nrows(); t++)
        {
            m_aif(t + 1) = ConcentrationFromSignal(aif_sig(t + 1), aif_sig0, aif_t1, aif_hct);
            //cerr << "t=" << t << "aif C=" << aif(t+1) << endl;
        }
    }

    // Case where we have a concentration signal already. No need to do any conversion.
    // Just copy the input AIF file
    else if (m_aif_type == "conc")
    {
        m_aif = read_ascii_matrix(args.GetString("aif-file"));
    }

    else
    {
        throw InvalidOptionValue("aif", m_aif_type, "Must be signal, conc or orton");
    }

    // Optional parameters. All of these can be inferred as well
    //m_vp = args.GetDoubleDefault("vp", 0);
    m_T10 = args.GetDoubleDefault("t10", 1);
    m_sig0 = args.GetDoubleDefault("sig0", 1000);
    m_delay = args.GetDoubleDefault("delay", 0);
    m_infer_vp = args.ReadBool("infer-vp");
    m_infer_t10 = args.ReadBool("infer-t10");
    m_infer_sig0 = args.ReadBool("infer-sig0");
    m_infer_delay = args.ReadBool("infer-delay");
    m_conv_method = args.GetStringDefault("convolution-method", "normal");
    m_use_log = args.ReadBool("use-log-params");
    //cerr << "Done Init" << endl;
}

vector<string> DCE_CTU_FwdModel::GetUsage() const
{
    vector<string> usage;

    usage.push_back("\nThis model implements the non-linear least square solution of the two compartment exchange model\n");
    return usage;
}

std::string DCE_CTU_FwdModel::GetDescription() const
{
    return "Non-linear least square solution of the two compartment exchange model";
}

void DCE_CTU_FwdModel::DumpParameters(const ColumnVector &vec, const string &indent) const
{
}

void DCE_CTU_FwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    names.push_back("Fp");
    names.push_back("PS");
    names.push_back("Vp");
    //names.push_back("Ve");
    if (m_infer_t10)
        names.push_back("T10");
    if (m_infer_sig0)
        names.push_back("sig0");
    if (m_infer_delay)
        names.push_back("delay");
}

double DCE_CTU_FwdModel::LogOrNot(double p) const
{
    if (m_use_log)
        return log(p);
    else
        return p;
}

void DCE_CTU_FwdModel::HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const
{
    //cerr << "priors" << endl;
    assert(prior.means.Nrows() == NumParams());

    SymmetricMatrix precs_prior = IdentityMatrix(NumParams()) * 1e-12;
    SymmetricMatrix precs_post = IdentityMatrix(NumParams()) * 1e-12;
    int p = 1;

    // Fp
    prior.means(p) = LogOrNot(initial_Fp);
    precs_prior(p, p) = 0.01;
    precs_post(p, p) = 0.1;
    p++;

    // PS
    prior.means(p) = LogOrNot(initial_PS);
    precs_prior(p, p) = 0.01;
    precs_post(p, p) = 0.1;
    p++;

    // Vp
    prior.means(p) = LogOrNot(initial_Vp);
    precs_prior(p, p) = 0.01;
    precs_post(p, p) = 1;
    p++;

    // We don't have Ve in this model
    // Ve
    //prior.means(p) = LogOrNot(0.02);
    //precs_prior(p, p) = 1e-20;
    //precs_post(p, p) = 1;
    //p++;

    // T10
    if (m_infer_t10)
    {
        prior.means(p) = LogOrNot(m_T10);
        precs_prior(p, p) = 1;
        precs_post(p, p) = 0.1;
        p++;
    }

    // sig0
    if (m_infer_sig0)
    {
        prior.means(p) = m_sig0;
        precs_prior(p, p) = 1e-12;
        precs_post(p, p) = 0.1;
        p++;
    }

    // delay
    if (m_infer_delay)
    {
        prior.means(p) = 0;
        precs_prior(p, p) = 1e-2;
        precs_post(p, p) = 1e-2;
        p++;
    }

    // Set precsions on priors
    prior.SetPrecisions(precs_prior);

    // Set initial posterior
    posterior = prior;
    posterior.SetPrecisions(precs_post);
    //cerr << "done priors" << endl;
}


ColumnVector DCE_CTU_FwdModel::aifshift(const ColumnVector &aif, const double delay) const
{
    // Shift a vector in time by interpolation (linear)
    // NB Makes assumptions where extrapolation is called for.

    // Number of whole time points of shift.
    int nshift = floor(delay / m_dt);
    // Fractional part of the shift
    double fshift = (delay / m_dt) - nshift;

    if (m_aif.Nrows() != data.Nrows())
        throw InvalidOptionValue("aif", stringify(m_aif.Nrows()), "Must have " + stringify(data.Nrows()) + " rows to match input data");
    ColumnVector aifnew(m_aif);
    int index;
    for (int i = 1; i <= data.Nrows(); i++)
    {
        index = i - nshift;
        if (index == 1)
        {
            //linear interpolation with zero as 'previous' time point. Only
            // possible if delay is > 0, so fshift > 0
            aifnew(i) = aif(1) * (1 - fshift);
        }
        else if (index < 1)
        {
            // Assume aif is zero before zeroth time point
            aifnew(i) = 0;
        }
        else if (index > data.Nrows())
        {
            // Beyond the final time point - assume aif takes the value of the final time point
            aifnew(i) = aif(data.Nrows());
        }
        else
        {
            //linear interpolation
            aifnew(i) = aif(index) + (aif(index - 1) - aif(index)) * fshift;
        }
    }
    return aifnew;
}

// This is the MR signal of spoiled gradient echo sequence
// Needs checking
double DCE_CTU_FwdModel::SignalFromConcentration(double C, double t10, double current_sig0) const
{
    double R10 = 1 / t10;
    double R1 = m_r1 * C + R10;
    double A = exp(-m_TR * R1);

    return current_sig0 * sin(m_FA) * (1 - A) / (1 - cos(m_FA) * A);
}

// Appendix A in https://www.sciencedirect.com/science/article/pii/S0730725X10001748
// This is converting the signal of AIF into concentration values
// Needs checking
double DCE_CTU_FwdModel::ConcentrationFromSignal(double s, double s0, double t10, double hct) const
{
    double e10 = exp(-m_TR / t10);
    double b = (1 - e10) / (1.0 - cos(m_FA) * e10);
    double sa = s / s0;
    double v = -log((1 - sa * b) / (1 - sa * b * cos(m_FA)));
    double r1 = v / m_TR;
    return ((r1 - 1 / t10) / m_r1) / (1 - hct);
}

// Implementation of Equation 4 of the paper
ColumnVector DCE_CTU_FwdModel::compute_concentration(double Fp, double PS, double Vp, const NEWMAT::ColumnVector &aif) const
{
    if(Fp == 0) {
        Fp = 0.00001;
    }
    if(PS == 0) {
        PS = 0.00001;
    }
    if(Vp == 0) {
        Vp = 0.00001;
    }

    double Tp = Vp / (Fp + PS);
    double E = PS / (PS + Fp);
    double Ktrans = E * Fp;

    // Result concentration
    ColumnVector current_concentration(data.Nrows());

    if(m_conv_method == "normal") {
        // Compute convolution using normal technique
        ColumnVector convolution_result(data.Nrows());

        // Bracket term in Equation 4
        ColumnVector bracket_term(data.Nrows());

        for (int t_index = 0; t_index < data.Nrows(); t_index++) {
            double current_t_value = t_index * m_dt - m_delay;
            bracket_term(t_index + 1) = Fp * exp(-current_t_value / Tp) + Ktrans * (1 - exp(-current_t_value / Tp));
            current_concentration(t_index + 1) = 0.0; // Initialize the concentration result vector. It seems that C++ likes this. :O
        }

        // Do convolution. 
        convolution_result = compute_convolution_normal(aif, bracket_term);
        // Compute concentration
        current_concentration = convolution_result;
    }

    // Need to work on this (you need to re-write the equation)
    if(m_conv_method == "iterative") {
        // Compute convolution using iterative technique 
        // We use the method in the appendix of this paper to compute the convolution
        //Some properties of convolution. X represents convolution operation
        //f X g = g X f
        //a * (f X g) = (a * f X g)
        //(f + h) X g = f X g + h X g
        //We need to rearrange the terms in Equation 4 using these properties
        //"""We need to rearrange Equation 4 of the model paper into the following way:"""
        //"""C(t) = (Fp - Ktrans) * Tp * (exp(-t/Tp))/Tp X Ca + Ktrans * Ca"""
        current_concentration = (Fp - Ktrans) * Tp * compute_convolution_iterative(aif, Tp) + Ktrans * aif;
    }

    return current_concentration;
}

// Compute convolution using normal method
// Assuming that the two input vectors have the same length
// https://stackoverflow.com/questions/24518989/how-to-perform-1-dimensional-valid-convolution
ColumnVector DCE_CTU_FwdModel::compute_convolution_normal(const NEWMAT::ColumnVector &term_1, const NEWMAT::ColumnVector &term_2) const
{
    // We only need the first Nrows() terms
    ColumnVector convolution_result(term_1.Nrows());

    for (int t_index = 0; t_index < term_1.Nrows(); t_index++) {
        convolution_result(t_index + 1) = 0;
        int jmn = (t_index >= term_1.Nrows() - 1)? t_index - (term_1.Nrows() - 1) : 0;
        int jmx = (t_index <  term_1.Nrows() - 1)? t_index                      : term_1.Nrows() - 1;
        for (int j = jmn; j <= jmx; j++) {
            convolution_result(t_index + 1) += (term_1(j + 1) * term_2(t_index - j + 1));
        }
    }
    return convolution_result;
}

// Implementation of appendix of the paper
ColumnVector DCE_CTU_FwdModel::compute_convolution_iterative(const NEWMAT::ColumnVector &aif, const double T_term) const
{
    int nhtpts = aif.Nrows();
    ColumnVector f_vector(nhtpts);

    if(T_term == 0) {
        f_vector = aif;
    }
    else {
        ColumnVector E0_vector(nhtpts - 1);
        ColumnVector E1_vector(nhtpts - 1);
        ColumnVector xi_vector(nhtpts - 1);
        ColumnVector aif_prime_vector(nhtpts - 1);

        for (int i = 2; i <= nhtpts; i++) {
            xi_vector(i - 1) = m_dt / T_term; // Because in the equation (A5) t_i+1 - t_i is always m_dt
            aif_prime_vector(i - 1) = (aif(i) - aif(i - 1)) / m_dt;
            E0_vector(i - 1) = 1 - exp(-xi_vector(i - 1));
            E1_vector(i - 1) = xi_vector(i - 1) - E0_vector(i - 1);
        }

        f_vector(1) = 0.0;

        for (int i = 2; i <= nhtpts; i++) {
            f_vector(i) = exp(-xi_vector(i - 1)) * f_vector(i - 1) + aif(i - 1) * E0_vector(i - 1) + aif_prime_vector(i - 1) * T_term * E1_vector(i - 1);
        }

    }

    return f_vector;
}


void DCE_CTU_FwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    //cerr << "evaluate" << endl;
    // Convert values if inferring logs, otherwise check not negative
    ColumnVector paramcpy = params;
    for (int i = 0; i < paramcpy.Nrows(); i++)
    {
        //cerr << "param " << i << "=" << paramcpy(i+1) << endl;
        if (m_use_log)
            paramcpy(i + 1) = exp(params(i + 1));
        else if (params(i + 1) < 0)
            paramcpy(i + 1) = 0;
    }

    // parameters that are inferred - extract and give sensible names
    //double Ktrans = paramcpy(1);
    //double Ve = paramcpy(2);
    //if (Ve == 0) cerr << "WARNING: Ve=0"<< endl;
    int p = 1;
    double Fp = paramcpy(p);
    ++p;
    double PS = paramcpy(p);
    ++p;
    double Vp = paramcpy(p);
    ++p;
    //double Ve = paramcpy(p); There is no Ve paramter in this model
    //double Vp = m_vp;
    double T10 = m_T10;
    double sig0 = m_sig0;
    double delay = m_delay;
    /*
    if (m_infer_vp)
    {
        Vp = paramcpy(p);
        p++;
    }
    */
    if (m_infer_t10)
    {
        T10 = paramcpy(p);
        p++;
    }
    if (m_infer_sig0)
    {
        sig0 = paramcpy(p);
        p++;
    }
    if (m_infer_delay)
    {
        delay = paramcpy(p);
        p++;
    }

    //ColumnVector C;
    ColumnVector concentration_tissue; // Tissue concentration results

    ColumnVector aif_current = m_aif;
    // Condition where we need to shift AIF
    if (m_infer_delay) {
        ColumnVector aif_current = aifshift(m_aif, delay);
    }

    concentration_tissue = compute_concentration(Fp, PS, Vp, aif_current);


    for (int i = 1; i <= concentration_tissue.Nrows(); i++)
    {
        if (isnan(concentration_tissue(i)) || isinf(concentration_tissue(i)))
        {
            cerr << "Warning NaN or inf in concentration_tissue" << endl;
            cerr << "result: " << concentration_tissue.t() << endl;
            cerr << "params: " << params.t() << endl;
            break;
        }
    }

    result.ReSize(data.Nrows());
    // Converts concentration back to MRI signal
    for (int i = 1; i <= data.Nrows(); i++)
    {   
        result(i) = SignalFromConcentration(concentration_tissue(i), T10, sig0);
    }

    for (int i = 1; i <= data.Nrows(); i++)
    {
        if (isnan(result(i)) || isinf(result(i)))
        {
            LOG << "Warning NaN or inf in result" << endl;
            LOG << "result: " << result.t() << endl;
            LOG << "params: " << params.t() << endl;
            LOG << "used params: " << paramcpy.t() << endl;

            result = 0.0;
            //getchar();
            break;

        }
    }
}

FwdModel *DCE_CTU_FwdModel::NewInstance()
{
    return new DCE_CTU_FwdModel();
}
