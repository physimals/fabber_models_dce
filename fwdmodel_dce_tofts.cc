/*  fwdmodel_dce.cc - Implements a convolution based model for DCE analysis

 Jesper Kallehauge, IBME

 Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce_tofts.h"

#include "newimage/newimageall.h"
#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "fabber_core/easylog.h"
#include "miscmaths/miscprob.h"

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DCEStdToftsFwdModel> DCEStdToftsFwdModel::registration("dce_tofts");

static OptionSpec OPTIONS[] = {
    { "delt", OPT_FLOAT, "Time resolution between volumes, in minutes", OPT_REQ, "" },
    { "fa", OPT_FLOAT, "Flip angle in degrees.", OPT_REQ, "" },
    { "tr", OPT_FLOAT, "Repetition time (TR) In seconds.", OPT_REQ, "" },
    { "r1", OPT_FLOAT, "Relaxivity of contrast agent, In s^-1 mM^-1.", OPT_REQ, "" },
    { "aif", OPT_STR, "Source of AIF function: orton=Orton (2008) population AIF, sig=User-supplied vascular signal, conc=User-supplied concentration curve", OPT_REQ, "none" },
    { "sig0", OPT_FLOAT, "Baseline signal. May be inferred.", OPT_NONREQ, "1000" },
    { "t10", OPT_FLOAT, "Baseline T1 value in seconds. May be inferred.", OPT_NONREQ, "1" },
    { "delay", OPT_FLOAT, "Delay time offset relative to AIF in minutes. May be inferred.", OPT_NONREQ, "0" },
    { "vp", OPT_FLOAT, "Fractional volume of blood plasma (Vp) in tissue", OPT_NONREQ, "0" },
    { "infer-vp", OPT_BOOL, "Infer the Vp parameter", OPT_NONREQ, "" },
    { "infer-delay", OPT_BOOL, "Infer the delay parameter", OPT_NONREQ, "" },
    { "infer-sig0", OPT_BOOL, "Infer baseline signal", OPT_NONREQ, "" },
    { "infer-T10", OPT_BOOL, "Infer T10 value", OPT_NONREQ, "" },
    { "use-log-params", OPT_BOOL, "Infer log values of parameters", OPT_NONREQ, "" },
    { "aif-file", OPT_FILE,
        "File containing single-column ASCII data defining the AIF. For aif=signal, this is the vascular signal curve. For aif=conc, it should be the blood plasma concentration curve",
        OPT_NONREQ, "none" },
    { "aif-hct", OPT_FLOAT, "Haematocrit value to use when converting an AIF signal to concentration. Used when aif=sig", OPT_NONREQ, "0.45"},
    { "aif-t1b", OPT_FLOAT, "Blood T1 value to use when converting an AIF signal to concentration. Used when aif=sig", OPT_NONREQ, "1.4"},
    { "aif-ab", OPT_FLOAT, "aB parameter for Orton AIF in mM. Used when aif=orton", OPT_NONREQ, "2.84"},
    { "aif-ag", OPT_FLOAT, "aG parameter for Orton AIF in min^-1. Used when aif=orton", OPT_NONREQ, "1.36"},
    { "aif-mub", OPT_FLOAT, "MuB parameter for Orton AIF in min^-1. Used when aif=orton", OPT_NONREQ, "22.8"},
    { "aif-mug", OPT_FLOAT, "MuG parameter for Orton AIF in min^-1. Used when aif=orton", OPT_NONREQ, "0.171"},
    { "" },
};

void DCEStdToftsFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string DCEStdToftsFwdModel::ModelVersion() const
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

void DCEStdToftsFwdModel::Initialize(FabberRunData &args)
{
    //cerr << "Initialize" << endl;
    // Mandatory parameters
    m_dt = args.GetDouble("delt", 0);
    m_FA = args.GetDouble("fa", 0, 90) * 3.1415926 / 180;
    m_TR = args.GetDouble("tr");
    m_r1 = args.GetDouble("r1");
    m_aif_type = args.GetString("aif");
    if (m_aif_type == "signal") {
        ColumnVector aif_sig = read_ascii_matrix(args.GetString("aif-file"));
        // Convert AIF to concentration curve
        float aif_t1 = args.GetDoubleDefault("aif-t1", 1.4);
        float aif_hct = args.GetDoubleDefault("aif-hct", 0.45);
        float aif_sig0 = aif_sig(1);
        m_aif = aif_sig;
        for (int t=0; t<m_aif.Nrows(); t++) {
            m_aif(t+1) = ConcentrationFromSignal(aif_sig(t+1), aif_sig0, aif_t1, aif_hct);
            //cerr << "t=" << t << "aif C=" << aif(t+1) << endl;
        }
    }
    else if (m_aif_type == "conc") {
        m_aif = read_ascii_matrix(args.GetString("aif-file"));
    }
    else if (m_aif_type == "orton") {
        m_ab = args.GetDoubleDefault("aif-ab", 2.84);
        m_ag = args.GetDoubleDefault("aif-ag", 1.36);
        m_mub = args.GetDoubleDefault("aif-mub", 22.8);
        m_mug = args.GetDoubleDefault("aif-mug", 0.171);
    }
    else {
        throw InvalidOptionValue("aif", m_aif_type, "Must be signal, conc or orton");
    }

    // Optional parameters. All of these can be inferred as well
    m_vp = args.GetDoubleDefault("vp", 0);
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

vector<string> DCEStdToftsFwdModel::GetUsage() const
{
    vector<string> usage;

    usage.push_back("\nThis model implements the standard and Extended Tofts one compartment model\n");
    return usage;
}

std::string DCEStdToftsFwdModel::GetDescription() const
{
    return "Standard and Extended Tofts one compartment model";
}

void DCEStdToftsFwdModel::DumpParameters(const ColumnVector &vec, const string &indent) const
{
}

void DCEStdToftsFwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    names.push_back("Ktrans");
    names.push_back("Ve");
    if (m_infer_vp) names.push_back("Vp");
    if (m_infer_t10) names.push_back("T10");
    if (m_infer_sig0) names.push_back("sig0");
    if (m_infer_delay) names.push_back("delay");
}

float DCEStdToftsFwdModel::LogOrNot(float p) const
{
    if (m_use_log) return log(p);
    else return p;
}

void DCEStdToftsFwdModel::HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const
{
    //cerr << "priors" << endl;
    assert(prior.means.Nrows() == NumParams());

    SymmetricMatrix precs_prior = IdentityMatrix(NumParams()) * 1e-12;
    SymmetricMatrix precs_post = IdentityMatrix(NumParams()) * 1e-12;
    int p = 1;

    // Ktrans
    prior.means(p) = LogOrNot(0.1);
    precs_prior(p,p) = 1e-20;
    precs_post(p,p) = 1;
    p++;
    
    // Ve
    prior.means(p) = LogOrNot(0.1);
    precs_prior(p,p) = 1e-20;
    precs_post(p,p) = 1;
    p++;

    // Vp
    if (m_infer_vp) {
        prior.means(p) = LogOrNot(0.01);
        precs_prior(p,p) = 100;
        precs_post(p,p) = 100;
        p++;
    }

    // T10
    if (m_infer_t10) {
        prior.means(p) = LogOrNot(m_T10);
        precs_prior(p,p) = 1e-12;
        precs_post(p,p) = 0.1;
        p++;
    }

    // sig0
    if (m_infer_sig0) {
        prior.means(p) = m_sig0;
        precs_prior(p,p) = 1e-12;
        precs_post(p,p) = 0.1;
        p++;
    }

    // delay
    if (m_infer_delay)
    {
        prior.means(p) = 0;
        precs_prior(p,p) = 1e-2;
        precs_post(p,p) = 1e-2;
        p++;
    }

    // Set precsions on priors
    prior.SetPrecisions(precs_prior);

    // Set initial posterior
    posterior = prior;
    posterior.SetPrecisions(precs_post);
    //cerr << "done priors" << endl;
}

ColumnVector DCEStdToftsFwdModel::GetConcentrationMeasuredAif(float delay, float Vp, float Ktrans, float Ve) const
{
    ColumnVector myaif = m_aif;
    if (m_infer_delay) ColumnVector myaif = aifshift(m_aif, delay);
    ColumnVector f(data.Nrows());

    float kep = Ktrans / Ve;
    for (int t = 0; t < data.Nrows(); t++)
    {
        float c = Vp * myaif(t+1);
        // Convolution integral using trapezium rule: I=INTEGRAL{aif(t)*exp(-kep(t-tau)*dt) }
        float I = (myaif(1)*exp(-kep*t*m_dt) + myaif(t+1))/2;
        //cerr << "t=" << t << ", c=" << c << endl;
        for (int tau=1; tau<t; tau++) {
            I += myaif(tau+1) * exp(-kep*(t-tau)*m_dt);
            //cerr << "t=" << t << ", tau=" << tau << ", c=" << c << endl;
        }
        I = I * Ktrans * m_dt;
        f(t+1) = c + I;
    }

    return f;
}

static float orton_f(float t, float a, float mub) 
{
    if (a == 0) return 0;

    float ret = ((1-exp(-a*t))/a) - (a*cos(mub*t) + mub*sin(mub*t) - a*exp(-a*t))/(a*a + mub*mub);
    return ret;
}

ColumnVector DCEStdToftsFwdModel::GetConcentrationOrton(float Vp, float Ktrans, float Ve) const
{
    ColumnVector f(data.Nrows());

    float kep = Ktrans / Ve;
    if (Ve == 0) {
        //cerr << "Ve=0 - returning 0" << endl;
        for (int tp = 0; tp < data.Nrows(); tp++) {
            f(tp+1) = 0;
        }
        return f;
    }

    for (int tp = 0; tp < data.Nrows(); tp++)
    {
        float c = 0;
        float t = tp*m_dt - m_delay;
        float tb = 2*3.14159265/m_mub;
        float Cp=0;
        if (t <= 0) {
            c = 0;
        }  
        else if (t <= tb) {
            Cp=m_ab*(1-cos(m_mub*t)) + m_ab*m_ag*orton_f(t, m_mug, m_mub);
            //LOG << "t<=tb=" << tb << ", Cp=" << Cp << endl;
            c = Vp * Cp;
            c += m_ab*m_ag*Ktrans * (orton_f(t, m_mug, m_mub) + (((kep-m_mug)/m_ag) - 1)*orton_f(t, kep, m_mub)) / (kep-m_mug);
            //cerr << "t=" << t << " : c=" << c << ", Cp=" << Cp << endl;
        }   
        else {
            Cp=m_ab*m_ag*orton_f(tb, m_mug, m_mub)*exp(-m_mug*(t-tb));
            //LOG << "t>tb=" << tb << ", Cp=" << Cp << endl;
            c = Vp * Cp;
            c += m_ab*m_ag*Ktrans * (orton_f(tb, m_mug, m_mub)*exp(-m_mug*(t-tb)) + (((kep-m_mug)/m_ag) - 1)*orton_f(tb, kep, m_mub)*exp(-kep*(t-tb))) / (kep-m_mug);
            //cerr << "t=" << t << " : c=" << c << ", Cp=" << Cp << endl;
        }     
        //cerr << " - t=" << tp << ", " << t << " : c=" << c << endl;
        if (isnan(c) || isinf(c)) {
            LOG << Ktrans << ", " << Ve << endl;
            LOG << "t=" << tp << " (" << t << ") tb=" << tb << ", Cp=" << Cp << ", Vp=" << Vp << ", kep=" << kep << ", mug=" << m_mug << ", ag=" << m_ag << " : c=" << c << endl;
            LOG << "term 1: " << m_ab*m_ag*Ktrans << endl;
            LOG << "term 2: " << orton_f(t, m_mug, m_mub) << endl;
            LOG << "term 3: " << (((kep-m_mug)/m_ag) - 1) << endl;
            LOG << "term 4: " << orton_f(t, kep, m_mub) << endl;
            LOG << "term 4.1: " << exp(-kep*t) << endl;
            LOG << "term 4.1: " << ((1-exp(-kep*t))/kep) << endl;
            LOG << "term 4.2: " << kep*cos(m_mub*t) << endl;
            LOG << "term 4.3: " << m_mub*sin(m_mub*t) << endl;
            LOG << "term 4.4: " << kep*exp(-kep*t) << endl;
            LOG << "term 4.5: " << (kep*kep + m_mub*m_mub) << endl;
        }

        f(tp+1) = c;
    }

    return f;
}

ColumnVector DCEStdToftsFwdModel::aifshift(const ColumnVector &aif, const float delay) const
{
    // Shift a vector in time by interpolation (linear)
    // NB Makes assumptions where extrapolation is called for.

    // Number of whole time points of shift.
    int nshift = floor(delay / m_dt);  
    // Fractional part of the shift
    float fshift = (delay / m_dt) - nshift; 
    
    if (m_aif.Nrows() != data.Nrows()) throw InvalidOptionValue("aif", stringify(m_aif.Nrows()), "Must have " + stringify(data.Nrows()) + " rows to match input data");
    ColumnVector aifnew(m_aif);
    int index;
    for (int i = 1; i <= data.Nrows(); i++)
    {
        index = i - nshift;
        if (index == 1)
        {
            //linear interpolation with zero as 'previous' time point. Only 
            // possible if delay is > 0, so fshift > 0
            aifnew(i) = aif(1) * (1-fshift);
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

float DCEStdToftsFwdModel::SignalFromConcentration(float C, float t10, float m0) const
{
    float R10 = 1/t10;
    //cerr << "R10=" << R10 << endl;
    //cerr << "Rg=" << m_r1 << endl;
    float R1 = m_r1 * C + R10;
    //cerr << "R1=" << R1 << endl;
    float A = exp(-m_TR * R1);
    //cerr << "A=" << A << endl;

    return m0 * sin(m_FA) * (1 - A) / (1 - cos(m_FA) * A);
}

float DCEStdToftsFwdModel::ConcentrationFromSignal(float s, float s0, float t10, float hct) const
{
    float e10 = exp(-m_TR / t10);
    float b = (1-e10) / (1.0-cos(m_FA) * e10);
    float sa = s/s0;
    float v = -log((1-sa*b)/(1-sa*b*cos(m_FA)));
    float r1 = v/m_TR;
	return ((r1 - 1/t10)/m_r1)/(1-hct);
}

void DCEStdToftsFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    //cerr << "evaluate" << endl;
    // Convert values if inferring logs, otherwise check not negative
    ColumnVector paramcpy = params;
    for (int i=0; i<paramcpy.Nrows(); i++) {
        //cerr << "param " << i << "=" << paramcpy(i+1) << endl;
        if (m_use_log) paramcpy(i+1) = exp(params(i+1));
        else if (params(i+1) < 0) paramcpy(i+1) = 0;
    }

    // parameters that are inferred - extract and give sensible names
    float Ktrans = paramcpy(1);
    float Ve = paramcpy(2);
    //if (Ve == 0) cerr << "WARNING: Ve=0"<< endl;
    int p = 2;
    float Vp = m_vp;
    float T10 = m_T10;
    float sig0 = m_sig0;
    float delay = m_delay;
    if (m_infer_vp) {
        Vp = paramcpy(p);
        p++;
    }
    if (m_infer_t10) {
        T10 = paramcpy(p);
        p++;
    }
    if (m_infer_sig0) {
        sig0 = paramcpy(p);
        p++;
    }
    if (m_infer_delay) {
        delay = paramcpy(p);
        p++;
    }
    
    //std::cerr << "TR=" << m_TR << ", FA=" << m_FA << ", r1=" << m_r1 << endl;
    //std::cerr << "Ktrans=" << Ktrans << ", Ve=" << Ve << ", T10=" << T10 << ", sig0=" << sig0 << endl;

    ColumnVector C;
    if (m_aif_type == "orton") {
        C = GetConcentrationOrton(Vp, Ktrans, Ve);
    }
    else {
        C = GetConcentrationMeasuredAif(delay, Vp, Ktrans, Ve);
    }

    //if (Ve == 0) {
    //    cerr << "Ve=0 - got C" << endl;
    //    cerr << C.t() << endl;
    //}

    /*for (int i = 1; i <= C.Nrows(); i++)
    {
        if (isnan(C(i)) || isinf(C(i)))
        {
            cerr << "Warning NaN or inf in C" << endl;
            cerr << "result: " << C.t() << endl;
            cerr << "params: " << params.t() << endl;

            result = 0.0;
            break;
        }
    }*/

    // Convert to the DCE signal
    float E10 = exp(-m_TR/T10);
    float m0 = sig0 * (1-cos(m_FA)*E10)/(sin(m_FA)*(1-E10));
    result.ReSize(data.Nrows());
    for (int i = 1; i <= data.Nrows(); i++)
    {
        result(i) = SignalFromConcentration(C(i), T10, m0);
        //std::cerr << "t=" << i-1 << ", C=" << C(i) << ", S=" << result(i) << endl;
     }
    
    for (int i = 1; i <= data.Nrows(); i++)
    {
        if (isnan(result(i)) || isinf(result(i)))
        {
            LOG << "Warning NaN or inf in result" << endl;
            LOG << "result: " << result.t() << endl;
            LOG << "params: " << params.t() << endl;

            result = 0.0;
            break;
        }
    }
    //cerr << "done evaluate" << endl;
}

FwdModel *DCEStdToftsFwdModel::NewInstance()
{
    return new DCEStdToftsFwdModel();
}
