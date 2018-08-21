/*  fwdmodel_dce_AATH_NLLS.cc - An Adiabatic Approximation to the Tissue Homogeneity Model for Water Exchange in the Brain: I. Theoretical Derivation

http://journals.sagepub.com/doi/full/10.1097/00004647-199812000-00011

 Moss Zhao - IBME, Oxford

 Copyright (C) 2018 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce_AATH_NLLS.h"

#include "fabber_core/easylog.h"
#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include <math.h>

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DCEAATHNLLSFwdModel> DCEAATHNLLSFwdModel::registration("dce_aath_nlls");

static OptionSpec OPTIONS[] = {
    { "delt", OPT_FLOAT, "Time resolution between volumes, in minutes", OPT_REQ, "" },
    { "fa", OPT_FLOAT, "Flip angle in degrees.", OPT_REQ, "" },
    { "tr", OPT_FLOAT, "Repetition time (TR) In seconds.", OPT_REQ, "" },
    { "r1", OPT_FLOAT, "Relaxivity of contrast agent, In s^-1 mM^-1.", OPT_REQ, "" },
    { "aif_type", OPT_STR, "Type of AIF must be aif_type=signal: User-supplied vascular MRI signal or aif_type=conc:User-supplied concentration curve", OPT_REQ, "none" },
    { "sig0", OPT_FLOAT, "Baseline signal. May be inferred.", OPT_NONREQ, "1000" },
    { "t10", OPT_FLOAT, "Baseline T1 value in seconds. May be inferred.", OPT_NONREQ, "1" },
    { "delay", OPT_FLOAT, "Delay time offset relative to AIF in minutes. May be inferred.", OPT_NONREQ, "0" },
    { "infer-delay", OPT_BOOL, "Infer the delay parameter", OPT_NONREQ, "" },
    { "infer-sig0", OPT_BOOL, "Infer baseline signal", OPT_NONREQ, "" },
    { "infer-T10", OPT_BOOL, "Infer T10 value", OPT_NONREQ, "" },
    { "use-log-params", OPT_BOOL, "Infer log values of parameters", OPT_NONREQ, "" },
    { "aif-file", OPT_FILE,
        "File containing single-column ASCII data defining the AIF. For aif=signal, this is the vascular signal curve. For aif=conc, it should be the blood plasma concentration curve",
        OPT_NONREQ, "none" },
    { "aif-hct", OPT_FLOAT, "Haematocrit value to use when converting an AIF signal to concentration. Used when aif=signal", OPT_NONREQ, "0.45" },
    { "aif-t1b", OPT_FLOAT, "Blood T1 value to use when converting an AIF signal to concentration. Used when aif=signal", OPT_NONREQ, "1.4" },
    { "" },
};

void DCEAATHNLLSFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string DCEAATHNLLSFwdModel::ModelVersion() const
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

void DCEAATHNLLSFwdModel::Initialize(FabberRunData &args)
{
    //cerr << "Initialize" << endl;
    // Mandatory parameters
    m_dt = args.GetDouble("delt", 0);
    m_FA = args.GetDouble("fa", 0, 90) * M_PI / 180; // convert from degree to radians
    m_TR = args.GetDouble("tr");
    m_r1 = args.GetDouble("r1");
    m_aif_type = args.GetString("aif_type");

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
    m_use_log = args.ReadBool("use-log-params");
    //cerr << "Done Init" << endl;
}

vector<string> DCEAATHNLLSFwdModel::GetUsage() const
{
    vector<string> usage;

    usage.push_back("\nThis model implements the adiabatic approximation to the tissue homogeneity model\n");
    return usage;
}

std::string DCEAATHNLLSFwdModel::GetDescription() const
{
    return "Adiabatic Approximation to the Tissue Homogeneity Model";
}

void DCEAATHNLLSFwdModel::DumpParameters(const ColumnVector &vec, const string &indent) const
{
}

void DCEAATHNLLSFwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    names.push_back("Fp");
    names.push_back("PS");
    names.push_back("Vp");
    names.push_back("Ve");
    if (m_infer_t10)
        names.push_back("T10");
    if (m_infer_sig0)
        names.push_back("sig0");
    if (m_infer_delay)
        names.push_back("delay");
}

double DCEAATHNLLSFwdModel::LogOrNot(double p) const
{
    if (m_use_log)
        return log(p);
    else
        return p;
}

void DCEAATHNLLSFwdModel::HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const
{
    //cerr << "priors" << endl;
    assert(prior.means.Nrows() == NumParams());

    SymmetricMatrix precs_prior = IdentityMatrix(NumParams()) * 1e-12;
    SymmetricMatrix precs_post = IdentityMatrix(NumParams()) * 1e-12;
    int p = 1;

    // Fp
    prior.means(p) = LogOrNot(0.02);
    precs_prior(p, p) = 1e-20;
    precs_post(p, p) = 1;
    p++;

    // PS
    prior.means(p) = LogOrNot(0.02);
    precs_prior(p, p) = 1e-20;
    precs_post(p, p) = 1;
    p++;

    // Vp
    prior.means(p) = LogOrNot(0.02);
    precs_prior(p, p) = 1e-20;
    precs_post(p, p) = 1;
    p++;

    // Ve
    prior.means(p) = LogOrNot(0.02);
    precs_prior(p, p) = 1e-20;
    precs_post(p, p) = 1;
    p++;

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


ColumnVector DCEAATHNLLSFwdModel::aifshift(const ColumnVector &aif, const double delay) const
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
double DCEAATHNLLSFwdModel::SignalFromConcentration(double C, double t10, double current_sig0) const
{
    double R10 = 1 / t10;
    double R1 = m_r1 * C + R10;
    double A = exp(-m_TR * R1);

    return current_sig0 * sin(m_FA) * (1 - A) / (1 - cos(m_FA) * A);
}

// Appendix A in https://www.sciencedirect.com/science/article/pii/S0730725X10001748
// This is converting the signal of AIF into concentration values
// Needs checking
double DCEAATHNLLSFwdModel::ConcentrationFromSignal(double s, double s0, double t10, double hct) const
{
    double e10 = exp(-m_TR / t10);
    double b = (1 - e10) / (1.0 - cos(m_FA) * e10);
    double sa = s / s0;
    double v = -log((1 - sa * b) / (1 - sa * b * cos(m_FA)));
    double r1 = v / m_TR;
    return ((r1 - 1 / t10) / m_r1) / (1 - hct);
}

// Equation 14 of the paper
ColumnVector DCEAATHNLLSFwdModel::compute_concentration(double Fp, double PS, double Vp, double Ve, NEWMAT::ColumnVector &aif) const
{
    double E = 1 - exp(-PS / Fp); // Equation 4
    double k_adb = E * Fp / Ve; // Equation 15
    double Vi = Vp;

    // Bracket term in Equation 7
    ColumnVector exp_term(data.Nrows());

    ColumnVector convolution_result(data.Nrows());
    ColumnVector current_concentration(data.Nrows());

    for (int t_index = 0; t_index < data.Nrows(); t_index++) {
        double current_t_value = t_index * m_dt - m_delay;
        exp_term(t_index + 1) = exp(-k_adb * current_t_value); // Equation 14
    }

    // Do convolution. 
    convolution_result = compute_convolution_normal(aif, exp_term);
    // Compute concentration
    current_concentration = Vi * aif + E * Fp * convolution_result; // Equation 14

    return current_concentration;
}

// Compute convolution using normal method
// Assuming that the two input vectors have the same length
// https://stackoverflow.com/questions/24518989/how-to-perform-1-dimensional-valid-convolution
ColumnVector DCEAATHNLLSFwdModel::compute_convolution_normal(const NEWMAT::ColumnVector &term_1, const NEWMAT::ColumnVector &term_2) const
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


void DCEAATHNLLSFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
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
    double Fp = paramcpy(1);
    double PS = paramcpy(2);
    double Vp = paramcpy(3);
    double Ve = paramcpy(4);
    //int p = 2;
    int p = 5;
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

    concentration_tissue = compute_concentration(Fp, PS, Vp, Ve, aif_current);

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

    cerr << "concentration: " << concentration_tissue.t() << endl;
    cerr << "result: " << result.t() << endl;
    cerr << "data: " << data.t() << endl;
    cerr << "params: " << params.t() << endl;
    getchar();

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

FwdModel *DCEAATHNLLSFwdModel::NewInstance()
{
    return new DCEAATHNLLSFwdModel();
}
