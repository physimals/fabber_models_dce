/*  fwdmodel_dce.cc - Implements a convolution based model for DCE analysis

 Jesper Kallehauge, IBME

 Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce_tofts.h"

#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>

#include <miscmaths/miscprob.h>
#include <newimage/newimageall.h>

#include <newmatio.h>

#include <iostream>
#include <stdexcept>

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DCEStdToftsFwdModel> DCEStdToftsFwdModel::registration("dce_tofts");

static OptionSpec OPTIONS[] = {
    { "delt", OPT_FLOAT, "Time resolution between volumes, in minutes", OPT_REQ, "" },
    { "fa", OPT_FLOAT, "Flip angle in degrees.", OPT_REQ, "" },
    { "tr", OPT_FLOAT, "Repetition time (TR) In seconds.", OPT_REQ, "" },
    { "r1", OPT_FLOAT, "Relaxivity of contrast agent, In s^-1 mM^-1.", OPT_REQ, "" },
    { "aif", OPT_STR, "Source of AIF function: orton=Orton (2008) population AIF, signal=User-supplied vascular signal, conc=User-supplied concentration curve", OPT_REQ, "none" },
    { "sig0", OPT_FLOAT, "Baseline signal. This value is ignored if sig0 is inferred.", OPT_NONREQ, "1" },
    { "t10", OPT_FLOAT, "Baseline T1 value in seconds. May be inferred.", OPT_NONREQ, "1" },
    { "delay", OPT_FLOAT, "Injection time (or delay time when using measured AIF) in minutes. May be inferred.", OPT_NONREQ, "0" },
    { "vp", OPT_FLOAT, "Fractional volume of blood plasma (Vp) in tissue", OPT_NONREQ, "0" },
    { "infer-vp", OPT_BOOL, "Infer the Vp parameter", OPT_NONREQ, "" },
    { "infer-delay", OPT_BOOL, "Infer the delay parameter", OPT_NONREQ, "" },
    { "infer-ve", OPT_BOOL, "Infer Ve rather than kep. Normally inferring kep is more numerically stable.", OPT_NONREQ, "" },
    { "infer-sig0", OPT_BOOL, "Infer baseline signal", OPT_NONREQ, "" },
    { "infer-t10", OPT_BOOL, "Infer t10 value", OPT_NONREQ, "" },
    { "auto-init-delay", OPT_BOOL, "Automatically initialize posterior value of delay parameter", OPT_NONREQ, "" },
    { "aif-file", OPT_FILE,
        "File containing single-column ASCII data defining the AIF. For aif=signal, this is the vascular signal curve. For aif=conc, it should be the blood plasma concentration curve",
        OPT_NONREQ, "none" },
    { "aif-hct", OPT_FLOAT, "Haematocrit value to use when converting an AIF signal to concentration. Used when aif=sig", OPT_NONREQ, "0.45" },
    { "aif-t1b", OPT_FLOAT, "Blood T1 value to use when converting an AIF signal to concentration. Used when aif=sig", OPT_NONREQ, "1.4" },
    { "aif-ab", OPT_FLOAT, "aB parameter for Orton AIF in mM. Used when aif=orton", OPT_NONREQ, "2.84" },
    { "aif-ag", OPT_FLOAT, "aG parameter for Orton AIF in min^-1. Used when aif=orton", OPT_NONREQ, "1.36" },
    { "aif-mub", OPT_FLOAT, "MuB parameter for Orton AIF in min^-1. Used when aif=orton", OPT_NONREQ, "22.8" },
    { "aif-mug", OPT_FLOAT, "MuG parameter for Orton AIF in min^-1. Used when aif=orton", OPT_NONREQ, "0.171" },
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

void DCEStdToftsFwdModel::Initialize(FabberRunData &rundata)
{
    // Mandatory parameters
    m_dt = rundata.GetDouble("delt", 0);
    m_fa = rundata.GetDouble("fa", 0, 90) * 3.1415926 / 180;
    m_tr = rundata.GetDouble("tr");
    m_r1 = rundata.GetDouble("r1");
    m_aif_type = rundata.GetString("aif");
    if (m_aif_type == "signal")
    {
        // Convert AIF signal to concentration curve
        ColumnVector aif_sig = read_ascii_matrix(rundata.GetString("aif-file"));
        double aif_t1 = rundata.GetDoubleDefault("aif-t1", 1.4);
        double aif_hct = rundata.GetDoubleDefault("aif-hct", 0.45);
        double aif_sig0 = aif_sig(1);
        m_aif = aif_sig;
        for (int t = 0; t < m_aif.Nrows(); t++)
        {
            m_aif(t + 1) = ConcentrationFromSignal(aif_sig(t + 1), aif_sig0, aif_t1, aif_hct);
        }
    }
    else if (m_aif_type == "conc")
    {
        m_aif = read_ascii_matrix(rundata.GetString("aif-file"));
    }
    else if (m_aif_type == "orton")
    {
        m_ab = rundata.GetDoubleDefault("aif-ab", 2.84);
        m_ag = rundata.GetDoubleDefault("aif-ag", 1.36);
        m_mub = rundata.GetDoubleDefault("aif-mub", 22.8);
        m_mug = rundata.GetDoubleDefault("aif-mug", 0.171);
    }
    else
    {
        throw InvalidOptionValue("aif", m_aif_type, "Must be signal, conc or orton");
    }

    // Additional model parameters. All of these can be inferred if required
    m_vp = rundata.GetDoubleDefault("vp", 0);
    m_t10 = rundata.GetDoubleDefault("t10", 1);
    m_sig0 = rundata.GetDoubleDefault("sig0", 1);
    m_delay = rundata.GetDoubleDefault("delay", 0);
    m_infer_vp = rundata.ReadBool("infer-vp");
    m_infer_t10 = rundata.ReadBool("infer-t10");
    m_infer_sig0 = rundata.ReadBool("infer-sig0");
    m_infer_delay = rundata.ReadBool("infer-delay");

    // Infer Ve rather than kep. kep is generally a bit more stable
    m_infer_ve = rundata.ReadBool("infer-ve");

    // Automatically initialise delay posterior
    m_auto_init_delay = rundata.ReadBool("auto-init-delay");
}

std::string DCEStdToftsFwdModel::GetDescription() const
{
    return "Standard and Extended Tofts one compartment model";
}

void DCEStdToftsFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    int p=0;
    params.push_back(Parameter(p++, "ktrans", DistParams(0.5, 1), DistParams(0.5, 1), PRIOR_NORMAL, TRANSFORM_ABS()));
    
    if (m_infer_ve) {
        params.push_back(Parameter(p++, "ve", DistParams(0.5, 1), DistParams(0.5, 1), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    }
    else {
        params.push_back(Parameter(p++, "kep", DistParams(0.5, 1), DistParams(1, 1), PRIOR_NORMAL, TRANSFORM_ABS()));
    }

    if (m_infer_sig0)
        params.push_back(Parameter(p++, "sig0", DistParams(1, 1e20), DistParams(1, 100), PRIOR_NORMAL, TRANSFORM_ABS()));
    if (m_infer_delay)
        params.push_back(Parameter(p++, "delay", DistParams(m_delay, 1), DistParams(m_delay, 1), PRIOR_NORMAL, TRANSFORM_ABS()));
    if (m_infer_vp)
        params.push_back(Parameter(p++, "vp", DistParams(0.05, 1), DistParams(0.05, 1), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    if (m_infer_t10)
        params.push_back(Parameter(p++, "t10", DistParams(m_t10, 1), DistParams(m_t10, 1), PRIOR_NORMAL, TRANSFORM_ABS()));
}

static int fit_step(const ColumnVector &data)
{
    // MSC method - fit a step function
    double best_ssq = -1;
    int step_pos = 0;
    for (int pos=1; pos<data.Nrows(); pos++) {
        double mean_left = data.Rows(1, pos).Sum() / pos;
        //cerr << "left=" << data_mean.Rows(1, ti).t() << " sum=" << data_mean.Rows(1, ti).Sum() << endl;
        double mean_right = data.Rows(pos+1, data.Nrows()).Sum() / (data.Nrows() - pos);
        //cerr << "right=" << data_mean.Rows(ti+1, data_mean.Nrows()).t() << " sum=" << data_mean.Rows(ti+1, data_mean.Nrows()).Sum() << endl;
        //cerr << "TI: " << ti << " (" << tis(ti) << ") mean left=" << mean_left << ", mean right=" << mean_right << endl;
        double ssq = 0;
        for (int t=1; t<=pos; t++) {
            ssq += (data(t) - mean_left)*(data(t) - mean_left);
        }
        for (int t=pos+1; t<=data.Nrows(); t++) {
            ssq += (data(t) - mean_right)*(data(t) - mean_right);
        }
        //cerr << "SSQ=" << ssq << " (best=" << best_ssq << ")" << endl;
        if ((ssq < best_ssq) || (best_ssq < 0)) {
            best_ssq = ssq;
            step_pos = pos;
        }
    }

    return step_pos-1;
}

void DCEStdToftsFwdModel::InitVoxelPosterior(MVNDist &posterior) const
{
    int sig0_idx = 3, delay_idx=3;
    if (m_infer_sig0) {
        posterior.means(sig0_idx) = data(1);
        delay_idx++;
    }
    if (m_infer_delay && m_auto_init_delay) {
        posterior.means(delay_idx) = m_dt * fit_step(data);
    }
}

ColumnVector DCEStdToftsFwdModel::GetConcentrationMeasuredAif(double delay, double vp, double ktrans, double kep) const
{
    ColumnVector myaif = m_aif;
    if (m_infer_delay)
        ColumnVector myaif = aifshift(m_aif, delay);
    ColumnVector f(data.Nrows());

    for (int t = 0; t < data.Nrows(); t++)
    {
        double c = vp * myaif(t + 1);
        // Convolution integral using trapezium rule: I=INTEGRAL{aif(t)*exp(-kep(t-tau)*dt) }
        double I = (myaif(1) * exp(-kep * t * m_dt) + myaif(t + 1)) / 2;
        for (int tau = 1; tau < t; tau++)
        {
            I += myaif(tau + 1) * exp(-kep * (t - tau) * m_dt);
        }
        I = I * ktrans * m_dt;
        f(t + 1) = c + I;
    }

    return f;
}

static double orton_f(double t, double a, double mub)
{
    if (a == 0)
        return 0;

    double ret = ((1 - exp(-a * t)) / a) - (a * cos(mub * t) + mub * sin(mub * t) - a * exp(-a * t)) / (a * a + mub * mub);
    return ret;
}

ColumnVector DCEStdToftsFwdModel::GetConcentrationOrton(double delay, double vp, double ktrans, double kep) const
{
    ColumnVector f(data.Nrows());

    for (int tp = 0; tp < data.Nrows(); tp++)
    {
        double c = 0;
        double t = tp * m_dt - delay;
        double tb = 2 * 3.14159265 / m_mub;
        double Cp = 0;
        if (t <= 0)
        {
            c = 0;
        }
        else if (t <= tb)
        {
            Cp = m_ab * (1 - cos(m_mub * t)) + m_ab * m_ag * orton_f(t, m_mug, m_mub);
            c = vp * Cp;
            c += m_ab * m_ag * ktrans * (orton_f(t, m_mug, m_mub) + (((kep - m_mug) / m_ag) - 1) * orton_f(t, kep, m_mub)) / (kep - m_mug);
        }
        else
        {
            Cp = m_ab * m_ag * orton_f(tb, m_mug, m_mub) * exp(-m_mug * (t - tb));
            c = vp * Cp;
            c += m_ab * m_ag * ktrans * (orton_f(tb, m_mug, m_mub) * exp(-m_mug * (t - tb)) + (((kep - m_mug) / m_ag) - 1) * orton_f(tb, kep, m_mub) * exp(-kep * (t - tb))) / (kep - m_mug);
        }
        
        if (isnan(c) || isinf(c))
        {
            LOG << ktrans << ", " << kep << endl;
            LOG << "t=" << tp << " (" << t << ") tb=" << tb << ", Cp=" << Cp << ", vp=" << vp << ", kep=" << kep << ", mug=" << m_mug << ", ag=" << m_ag << " : c=" << c << endl;
            LOG << "term 1: " << m_ab * m_ag * ktrans << endl;
            LOG << "term 2: " << orton_f(t, m_mug, m_mub) << endl;
            LOG << "term 3: " << (((kep - m_mug) / m_ag) - 1) << endl;
            LOG << "term 4: " << orton_f(t, kep, m_mub) << endl;
            LOG << "term 4.1: " << exp(-kep * t) << endl;
            LOG << "term 4.1: " << ((1 - exp(-kep * t)) / kep) << endl;
            LOG << "term 4.2: " << kep * cos(m_mub * t) << endl;
            LOG << "term 4.3: " << m_mub * sin(m_mub * t) << endl;
            LOG << "term 4.4: " << kep * exp(-kep * t) << endl;
            LOG << "term 4.5: " << (kep * kep + m_mub * m_mub) << endl;
        }

        f(tp + 1) = c;
    }

    return f;
}

ColumnVector DCEStdToftsFwdModel::aifshift(const ColumnVector &aif, const double delay) const
{
    // Shift a vector in time by interpolation (linear)
    // NB Makes assumptions where extrapolation is called for.

    // Number of whole time points of shift.
    int nshift = int(floor(delay / m_dt));
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

double DCEStdToftsFwdModel::SignalFromConcentration(double C, double t10, double m0) const
{
    double R10 = 1 / t10;
    double R1 = m_r1 * C + R10;
    double A = exp(-m_tr * R1);

    return m0 * sin(m_fa) * (1 - A) / (1 - cos(m_fa) * A);
}

double DCEStdToftsFwdModel::ConcentrationFromSignal(double s, double s0, double t10, double hct) const
{
    double e10 = exp(-m_tr / t10);
    double b = (1 - e10) / (1.0 - cos(m_fa) * e10);
    double sa = s / s0;
    double v = -log((1 - sa * b) / (1 - sa * b * cos(m_fa)));
    double r1 = v / m_tr;
    return ((r1 - 1 / t10) / m_r1) / (1 - hct);
}

void DCEStdToftsFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Parameters that are inferred - extract and give sensible names
    double ktrans = params(1);
    double kep;
    if (m_infer_ve) {
        kep = ktrans / params(2);
    }
    else {
        kep = params(2);        
    }

    // Optional parameters to infer
    double sig0 = m_sig0;
    double delay = m_delay;
    double vp = m_vp;
    double t10 = m_t10;
    int p = 3;
    if (m_infer_sig0)
    {
        sig0 = params(p++);
    }
    if (m_infer_delay)
    {
        delay = params(p++);
    }
    if (m_infer_vp)
    {
        vp = params(p++);
    }
    if (m_infer_t10)
    {
        t10 = params(p++);
    }

    ColumnVector C;
    if (m_aif_type == "orton")
    {
        C = GetConcentrationOrton(delay, vp, ktrans, kep);
    }
    else
    {
        C = GetConcentrationMeasuredAif(delay, vp, ktrans, kep);
    }

    // Convert concentration curve to DCE signal
    double E10 = exp(-m_tr / t10);
    double m0 = sig0 * (1 - cos(m_fa) * E10) / (sin(m_fa) * (1 - E10));
    result.ReSize(data.Nrows());
    for (int i = 1; i <= data.Nrows(); i++)
    {
        result(i) = SignalFromConcentration(C(i), t10, m0);
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
}

FwdModel *DCEStdToftsFwdModel::NewInstance()
{
    return new DCEStdToftsFwdModel();
}
