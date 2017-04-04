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
    { "fa", OPT_FLOAT, "Flip angle in degrees. Required if Acq_tech is not none", OPT_REQ, "" },
    { "tr", OPT_FLOAT, "Repetition time (TR) In seconds. Required if Acq_tech is not none", OPT_REQ, "" },
    { "r1", OPT_FLOAT, "Relaxivity of contrast agent, In s^-1 mM^-1. Required if Acq_tech is not none", OPT_REQ, "" },
    { "aif", OPT_FILE,
        "Fixed AIF arterial signal as single-column ASCII data. Must be concentration curve, not signal curve",
        OPT_REQ, "none" },
    { "vp", OPT_FLOAT, "Fractional volume of blood plasma in tissue", OPT_NONREQ, "0" },
    { "delay", OPT_FLOAT, "Estimated injection time in minutes", OPT_NONREQ, "0" },
    { "t10", OPT_FLOAT, "Baseline T1 value in seconds. May be inferred if desired", OPT_NONREQ, "" },
    { "sig0", OPT_FLOAT, "Baseline signal. May be inferred if desired", OPT_NONREQ, "" },
    { "infer-vp", OPT_BOOL, "Infer the Vp parameter", OPT_NONREQ, "" },
    { "infer-delay", OPT_BOOL, "Infer the delay parameter", OPT_NONREQ, "" },
    { "infer-sig0", OPT_BOOL, "Infer baseline signal", OPT_NONREQ, "" },
    { "infer-T10", OPT_BOOL, "Infer T10 value", OPT_NONREQ, "" },
    { "use-log-params", OPT_BOOL, "Infer log values of parameters", OPT_NONREQ, "" },
    { "aif-is-signal", OPT_BOOL, "Interpret AIF as signal rather than concentration curve", OPT_NONREQ, "" },
    { "aif-hct", OPT_FLOAT, "Haemocrit value to use when converting AIF signal to concentration", OPT_NONREQ, "0.45"},
    { "aif-t1", OPT_FLOAT, "T1 value to use when converting AIF signal to concentration", OPT_NONREQ, "1.4"},
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

void DCEStdToftsFwdModel::Initialize(ArgsType &args)
{
    //cerr << "Initialize" << endl;
    // Mandatory parameters
    m_dt = convertTo<double>(args.Read("delt"));
    m_FA = convertTo<double>(args.GetString("fa")) * 3.1415926 / 180;;
    m_TR = convertTo<double>(args.GetString("tr"));
    m_r1 = convertTo<double>(args.GetString("r1"));
    string artfile = args.GetString("aif");
    aif = read_ascii_matrix(artfile);
    
    // Optional parameters. All of these can be inferred as well
    m_vp = convertTo<double>(args.GetStringDefault("vp", "0"));
    m_T10 = convertTo<double>(args.GetStringDefault("t10", "1"));
    m_sig0 = convertTo<double>(args.GetStringDefault("sig0", "1000"));
    m_delay = convertTo<double>(args.GetStringDefault("delay", "0"));
    m_infer_vp = args.ReadBool("infer-vp");
    m_infer_t10 = args.ReadBool("infer-t10");
    m_infer_sig0 = args.ReadBool("infer-sig0");
    m_infer_delay = args.ReadBool("infer-delay");

    m_use_log = args.ReadBool("use-log-params");
    
    ColumnVector aif_orig(aif);
    if (args.GetBool("aif-is-signal")) {
        // Convert AIF to concentration curve
        float aif_t1 = convertTo<double>(args.GetStringDefault("aif-t1", "1.4"));
        float aif_hct = convertTo<double>(args.GetStringDefault("aif-hct", "0.45"));
        float aif_sig0 = aif_orig(1);
        for (int t=0; t<aif.Nrows(); t++) {
            aif(t+1) = ConcentrationFromSignal(aif_orig(t+1), aif_sig0, aif_t1, aif_hct);
            //cerr << "t=" << t << "aif C=" << aif(t+1) << endl;
        }
    }
    //cerr << "Done Init" << endl;
}

vector<string> DCEStdToftsFwdModel::GetUsage() const
{
    vector<string> usage;

    usage.push_back("\nThis model is the standard Tofts one compartment model\n");
    return usage;
}

std::string DCEStdToftsFwdModel::GetDescription() const
{
    return "Standard Tofts one compartment model";
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
    cerr << "priors" << endl;
    assert(prior.means.Nrows() == NumParams());

    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;
    int p = 1;

    // Ktrans
    prior.means(p) = LogOrNot(0.1);
    precisions(p,p) = 1e-20;
    p++;
    
    // Ve
    prior.means(p) = LogOrNot(0.1);
    precisions(p,p) = 1e-20;
    p++;

    // Vp
    if (m_infer_vp) {
        prior.means(p) = LogOrNot(0.01);
        precisions(p,p) = 100;
        p++;
    }

    // T10
    if (m_infer_t10) {
        prior.means(p) = LogOrNot(m_T10);
        precisions(p,p) = 1e-12;
        p++;
    }

    // sig0
    if (m_infer_sig0) {
        prior.means(p) = m_sig0;
        precisions(p,p) = 1e-12;
        p++;
    }

    // delay
    if (m_infer_delay)
    {
        prior.means(p) = 0;
        precisions(p,p) = 1e-2;
        p++;
    }

    // Set precsions on priors
    prior.SetPrecisions(precisions);

    // Set initial posterior
    posterior = prior;
    posterior.SetPrecisions(precisions);
    cerr << "done priors" << endl;
}

ColumnVector DCEStdToftsFwdModel::convolution(const ColumnVector &myaif, float Vp, float Ktrans, float Ve) const
{
    int nt = myaif.Nrows();
    ColumnVector f(nt);

    float kep = Ktrans / Ve;
    for (int t = 0; t < nt; t++)
    {
        float c = Vp * myaif(t+1);
        // Convolution integral using trapezium rule
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

ColumnVector DCEStdToftsFwdModel::aifshift(const ColumnVector &aif, const float delay) const
{
    // Shift a vector in time by interpolation (linear)
    // NB Makes assumptions where extrapolation is called for.

    // Number of whole time points of shift.
    int nshift = floor(delay / m_dt);  
    // Fractional part of the shift
    float fshift = (delay / m_dt) - nshift; 
    
    ColumnVector aifnew(aif);
    int index;
    int nt = aif.Nrows();
    for (int i = 1; i <= nt; i++)
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
        else if (index > nt)
        {
            // Beyond the final time point - assume aif takes the value of the final time point
            aifnew(i) = aif(nt);
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
    // ensure that values are reasonable
    // negative check
    ColumnVector paramcpy = params;
    if (m_use_log) {
        for (int i=0; i<paramcpy.Nrows(); i++) {
            paramcpy(i+1) = exp(params(i+1));
        }
    }
    else {
        for (int i = 1; i <= NumParams(); i++)
        {
            if (params(i) < 0)
            {
                paramcpy(i) = 0;
            }
        }
    }

   //for (int i=0; i<paramcpy.Nrows(); i++) {
   //      cerr << "param " << i << "=" << paramcpy(i+1) << endl;
   // }

    ColumnVector myaif(aif);
    int nt = aif.Nrows();
    assert(nt == result.Nrows());

    // parameters that are inferred - extract and give sensible names
    float Ktrans = paramcpy(1);
    float Ve = paramcpy(2);
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
        myaif = aifshift(aif, delay);
        p++;
    }
    
    //std::cerr << "TR=" << m_TR << ", FA=" << m_FA << ", r1=" << m_r1 << endl;
    //std::cerr << "Ktrans=" << Ktrans << ", Ve=" << Ve << ", T10=" << T10 << ", s0=" << s0 << ", Kep=" << Kep << ", Vp=" << Vp << endl;

    ColumnVector C = convolution(myaif, Vp, Ktrans, Ve);

    //convert to the DCE signal
    
    float E10 = exp(-m_TR/T10);
    float m0 = sig0 * (1-cos(m_FA)*E10)/(sin(m_FA)*(1-E10));
    result.ReSize(nt);
    for (int i = 1; i <= nt; i++)
    {
        result(i) = SignalFromConcentration(C(i), T10, m0);
        //std::cerr << "t=" << i-1 << ", C=" << C(i) << ", S=" << result(i) << endl;
     }
    
    for (int i = 1; i <= nt; i++)
    {
        if (isnan(result(i)) || isinf(result(i)))
        {
            LOG << "Warning NaN of inf in result" << endl;
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
