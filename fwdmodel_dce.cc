#include "fwdmodel_dce.h"

#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>

#include <miscmaths/miscprob.h>
#include <newimage/newimageall.h>

#include <newmatio.h>

#include <iostream>
#include <stdexcept>

using namespace NEWMAT;

static OptionSpec OPTIONS[] = {
    { "delt", OPT_FLOAT, "Time resolution between volumes, in minutes", OPT_REQ, "" },
    { "fa", OPT_FLOAT, "Flip angle in degrees.", OPT_REQ, "" },
    { "tr", OPT_FLOAT, "Repetition time (TR) In seconds.", OPT_REQ, "" },
    { "r1", OPT_FLOAT, "Relaxivity of contrast agent, In s^-1 mM^-1.", OPT_REQ, "" },

    { "t10", OPT_FLOAT, "Baseline T1 value in seconds. May be inferred.", OPT_NONREQ, "1" },
    { "sig0", OPT_FLOAT, "Baseline signal. This value is ignored if sig0 is inferred.", OPT_NONREQ, "1" },
    { "delay", OPT_FLOAT, "Injection time (or delay time when using measured AIF) in minutes. May be inferred.", OPT_NONREQ, "0" },
    { "infer-t10", OPT_BOOL, "Infer t10 value", OPT_NONREQ, "" },
    { "infer-sig0", OPT_BOOL, "Infer baseline signal", OPT_NONREQ, "" },
    { "infer-delay", OPT_BOOL, "Infer the delay parameter", OPT_NONREQ, "" },

    { "aif", OPT_STR, "Source of AIF function: orton=Orton (2008) population AIF, parker=Parker (2006) population AIF, signal=User-supplied vascular signal, conc=User-supplied concentration curve", OPT_REQ, "none"},
    { "aif-data", OPT_MATRIX,
        "File containing single-column ASCII data defining the AIF. For aif=signal, this is the vascular signal curve. For aif=conc, it should be the blood plasma concentration curve",
        OPT_NONREQ, "none" },
    { "aif-hct", OPT_FLOAT, "Haematocrit value to use when converting an AIF signal to concentration. Used when aif=signal", OPT_NONREQ, "0.45" },
    { "aif-t1b", OPT_FLOAT, "Blood T1 value to use when converting an AIF signal to concentration. Used when aif=signal", OPT_NONREQ, "1.4" },
    { "aif-ab", OPT_FLOAT, "aB parameter for Orton AIF in mM. Used when aif=orton", OPT_NONREQ, "2.84" },
    { "aif-ag", OPT_FLOAT, "aG parameter for Orton AIF in min^-1. Used when aif=orton", OPT_NONREQ, "1.36" },
    { "aif-mub", OPT_FLOAT, "MuB parameter for Orton AIF in min^-1. Used when aif=orton", OPT_NONREQ, "22.8" },
    { "aif-mug", OPT_FLOAT, "MuG parameter for Orton AIF in min^-1. Used when aif=orton", OPT_NONREQ, "0.171" },

    { "auto-init-delay", OPT_BOOL, "Automatically initialize posterior value of delay parameter", OPT_NONREQ, "" },
    { "" },
};

void DCEFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string DCEFwdModel::ModelVersion() const
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

void DCEFwdModel::Initialize(FabberRunData &rundata)
{
    // Mandatory parameters
    m_dt = rundata.GetDouble("delt", 0);
    m_fa = rundata.GetDouble("fa", 0, 90) * 3.1415926 / 180;
    m_tr = rundata.GetDouble("tr");
    m_r1 = rundata.GetDouble("r1");

    // AIF
    m_aif_type = rundata.GetString("aif");
    if (m_aif_type == "signal")
    {
        // Convert AIF signal to concentration curve
        ColumnVector aif_sig = read_ascii_matrix(rundata.GetString("aif-data"));
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
        m_aif = read_ascii_matrix(rundata.GetString("aif-data"));
    }
    else if (m_aif_type == "orton")
    {
        m_ab = rundata.GetDoubleDefault("aif-ab", 2.84);
        m_ag = rundata.GetDoubleDefault("aif-ag", 1.36);
        m_mub = rundata.GetDoubleDefault("aif-mub", 22.8);
        m_mug = rundata.GetDoubleDefault("aif-mug", 0.171);
    }
    else if (m_aif_type == "parker")
    {
        // FIXME read in Parker AIF parameters rather than hardcoding?
    }
    else
    {
        throw InvalidOptionValue("aif", m_aif_type, "Must be signal, conc, parker or orton");
    }

    // Standard DCE model parameters. All of these can be inferred if required
    m_t10 = rundata.GetDoubleDefault("t10", 1);
    m_sig0 = rundata.GetDoubleDefault("sig0", 1);
    m_delay = rundata.GetDoubleDefault("delay", 0);
    m_infer_t10 = rundata.ReadBool("infer-t10");
    m_infer_sig0 = rundata.ReadBool("infer-sig0");
    m_infer_delay = rundata.ReadBool("infer-delay");

    // Automatically initialise delay posterior
    m_auto_init_delay = rundata.ReadBool("auto-init-delay");
}

void DCEFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    // Add standard DCE parameters to whatever's there
    m_sig0_idx = params.size();
    unsigned int p = m_sig0_idx;
    if (m_infer_sig0)
        params.push_back(Parameter(p++, "sig0", DistParams(m_sig0, 1e8), DistParams(m_sig0, 100), PRIOR_NORMAL, TRANSFORM_ABS()));
    if (m_infer_delay)
        params.push_back(Parameter(p++, "delay", DistParams(m_delay, 100), DistParams(m_delay, 1), PRIOR_NORMAL, TRANSFORM_ABS()));
    if (m_infer_t10)
        params.push_back(Parameter(p++, "t10", DistParams(m_t10, 1), DistParams(m_t10, 1), PRIOR_NORMAL, TRANSFORM_ABS()));
}

static int fit_step(const ColumnVector &data)
{
    // MSC method for initializing delay parameter.
    // We fit a step function to the DCE signal to estimate the delay. Otherwise the
    // delay can be problematic to infer.
    double best_ssq = -1;
    int step_pos = 0;
    for (int pos=1; pos<data.Nrows(); pos++) {
        double mean_left = data.Rows(1, pos).Sum() / pos;
        double mean_right = data.Rows(pos+1, data.Nrows()).Sum() / (data.Nrows() - pos);
        double ssq = 0;
        for (int t=1; t<=pos; t++) {
            ssq += (data(t) - mean_left)*(data(t) - mean_left);
        }
        for (int t=pos+1; t<=data.Nrows(); t++) {
            ssq += (data(t) - mean_right)*(data(t) - mean_right);
        }
        if ((ssq < best_ssq) || (best_ssq < 0)) {
            best_ssq = ssq;
            step_pos = pos;
        }
    }

    return step_pos-1;
}

void DCEFwdModel::InitVoxelPosterior(MVNDist &posterior) const
{
    int delay_idx = m_sig0_idx;
    if (m_infer_sig0) {
        // Note that sig0 is the fully relaxed signal - not the
        // raw MRI signal. So to estimate it from the data we 
        // need this correction.
        double sig0_data = data(1);
        double E10 = exp(-m_tr / m_t10);
        double sig0_relaxed = sig0_data * (1 - cos(m_fa) * E10) / (sin(m_fa) * (1 - E10));
        posterior.means(m_sig0_idx+1) = sig0_relaxed;
        delay_idx++;
    }
    if (m_infer_delay && m_auto_init_delay) {
        posterior.means(delay_idx+1) = m_dt * fit_step(data);
    }
}

/**
 * Get the AIF function taking into account the delay parameter
 *
 * Uses linear interpolation and makes assumptions where extrapolation is called for.
 */
ColumnVector DCEFwdModel::aifshift(const ColumnVector &aif, const double delay) const
{
    // Don't bother if there is no delay
    if (delay == 0) 
    {
        return aif;
    }

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
            // linear interpolation with zero as 'previous' time point. Only
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
double DCEFwdModel::SignalFromConcentration(double C, double t10, double current_sig0) const
{
    double R10 = 1 / t10;
    double R1 = m_r1 * C + R10;
    double A = exp(-m_tr * R1);

    return current_sig0 * sin(m_fa) * (1 - A) / (1 - cos(m_fa) * A);
}

// Appendix A in https://www.sciencedirect.com/science/article/pii/S0730725X10001748
// This is converting the signal of AIF into concentration values
// Needs checking
double DCEFwdModel::ConcentrationFromSignal(double s, double s0, double t10, double hct) const
{
    double e10 = exp(-m_tr / t10);
    double b = (1 - e10) / (1.0 - cos(m_fa) * e10);
    double sa = s / s0;
    double v = -log((1 - sa * b) / (1 - sa * b * cos(m_fa)));
    double r1 = v / m_tr;
    return ((r1 - 1 / t10) / m_r1) / (1 - hct);
}

// f(t, a) from Orton 2008 (Eq A.2)
double DCEFwdModel::OrtonF(double t, double a) const
{
    if (a == 0)
        return 0;

    double ret = ((1 - exp(-a * t)) / a) - (a * cos(m_mub * t) + m_mub * sin(m_mub * t) - a * exp(-a * t)) / (a * a + m_mub * m_mub);
    return ret;
}

// Orton 20008 (Eq 4) - AIF function in concentration
double DCEFwdModel::OrtonAIF(double t) const
{
    double tb = 2*3.14159265 / m_mub;
    if (t <= 0) 
    {
        return 0;
    }
    else if ( t <= tb) 
    {
        return m_ab * (1-cos(m_mub*t)) + m_ab*m_ag*OrtonF(t, m_mug);
    }
    else 
    {
        return m_ab*m_ag*OrtonF(tb, m_mug)*exp(-m_mug*(t-tb));
    }
}

// Parker et al 2006: Equation 1
double DCEFwdModel::ParkerAIF(double t) const
{
    double A1 = 0.809;  // mmol * min
    double T1 = 0.17046;  // min
    double sigma1 = 0.0563;  // min
    
    double A2 = 0.330; // mmol * min
    double T2 = 0.365;  // min
    double sigma2 = 0.132; // min
    
    double alpha = 1.050; // mmol
    double beta = 0.1685; // min-1
    double s = 38.078; // min-1
    double tau = 0.483; // min

    // This is the summation terms. Two Gaussians
    double sum_term_1 = A1 * exp(-(t - T1)*(t - T1) / (2*sigma1*sigma1)) / (sigma1*sqrt(2*M_PI));
    double sum_term_2 = A2 * exp(-(t - T2)*(t - T2) / (2*sigma2*sigma2)) / (sigma2*sqrt(2*M_PI));

    // This is the exponential modulated with a sigmoid function term of Equation 1
    double sigmoid_term = alpha * exp(-beta * t) / (1 + exp(-s * (t - tau)));

    // Put them together
    return sum_term_1 + sum_term_2 + sigmoid_term; // unit of result should be mmol
}

double DCEFwdModel::AIF(double t) const
{
    if (t <= 0) 
    {
        return 0;
    }
    else if (m_aif_type == "orton") 
    {
        return OrtonAIF(t);
    }
    else if (m_aif_type == "parker") 
    {
        return ParkerAIF(t);
    }
    else 
    {
        // Linearly interpolate provided AIF. Assume AIF is 
        // constant beyond the last point
        double t_idx = t / m_dt;
        int t_idx_0 = int(floor(t_idx));
        double t_idx_frac = t_idx - double(t_idx_0);
        if (t_idx_0 < m_aif.Nrows()-1)
        {
            return t_idx_frac*m_aif(t_idx_0+1) + (1-t_idx_frac)*m_aif(t_idx_0+2);
        }
        else 
        {
            return m_aif(m_aif.Nrows());
        }
    }   
}