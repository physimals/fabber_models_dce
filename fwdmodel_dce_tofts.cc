/**
 * fwdmodel_dce_tofts.cc
 *
 * Standard and Extended Tofts model for DCE analysis
 *
 * Martin Craig, IBME
 *
 * Copyright (C) 2008 University of Oxford
 */

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

std::string DCEStdToftsFwdModel::GetDescription() const
{
    return "Standard and Extended Tofts one compartment model";
}

static OptionSpec OPTIONS[] = {
    { "ktrans", OPT_FLOAT, "Initial value of ktrans parameter", OPT_NONREQ, "0.5" },
    { "ve", OPT_FLOAT, "Initial value of ve parameter", OPT_NONREQ, "0.5" },
    { "kep", OPT_FLOAT, "Initial value of kep parameter if using --infer-kep", OPT_NONREQ, "1" },
    { "vp", OPT_FLOAT, "Fractional volume of blood plasma (Vp) in tissue", OPT_NONREQ, "0" },
    { "infer-vp", OPT_BOOL, "Infer the Vp parameter", OPT_NONREQ, "" },
    { "infer-kep", OPT_BOOL, "Infer kep rather than Ve. Sometimes inferring kep is more numerically stable.", OPT_NONREQ, "" },
    { "force-conv", OPT_BOOL, "Force numerical convolution - don't use analytic solution for Orton AIF", OPT_NONREQ, "" },
    { "" },
};

void DCEStdToftsFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    DCEFwdModel::GetOptions(opts);
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

void DCEStdToftsFwdModel::Initialize(FabberRunData &rundata)
{
    DCEFwdModel::Initialize(rundata);

    m_ktrans = rundata.GetDoubleDefault("ktrans", 0.5);
    m_infer_kep = rundata.ReadBool("infer-kep");
    if (m_infer_kep)
    {
        m_ve = m_ktrans / rundata.GetDoubleDefault("kep", 1);
    }
    else
    {
        m_ve = rundata.GetDoubleDefault("ve", 0.5);
    }

    m_infer_vp = rundata.ReadBool("infer-vp");
    m_vp = rundata.GetDoubleDefault("vp", m_infer_vp ? 0.05 : 0);

    m_force_conv = rundata.ReadBool("force-conv");
}

void DCEStdToftsFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    // Model parameters
    int p=0;
    params.push_back(Parameter(p++, "ktrans", DistParams(m_ktrans, 1e5), DistParams(m_ktrans, 100), PRIOR_NORMAL, TRANSFORM_LOG()));

    if (m_infer_kep)
    {
        double kep = m_ktrans / m_ve;
        params.push_back(Parameter(p++, "kep", DistParams(kep, 1e5), DistParams(kep, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
    }
    else
    {
        params.push_back(Parameter(p++, "ve", DistParams(m_ve, 1e5), DistParams(m_ve, 1), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    }

    if (m_infer_vp)
    {
        params.push_back(Parameter(p++, "vp", DistParams(m_vp, 1), DistParams(m_vp, 1), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    }

    // Standard DCE parameters
    DCEFwdModel::GetParameterDefaults(params);
}

ColumnVector DCEStdToftsFwdModel::GetConcentrationMeasuredAif(double delay, double vp, double ktrans, double kep) const
{
    ColumnVector f(data.Nrows());

    for (int t_idx = 0; t_idx < data.Nrows(); t_idx++)
    {
        double t = t_idx * m_dt - delay;
        if (t <= 0)
        {
            f(t_idx + 1) = 0;
        }
        else
        {
            // Convolution integral: I=INTEGRAL{aif(t)*exp(-kep(t-tau)*dt) }
            double I = 0;
            for (int tau_idx = 0; tau_idx < t_idx; tau_idx++)
            {
                double tau = tau_idx * m_dt - delay;
                if (tau > 0) {
                    I += AIF(tau) * exp(-kep * (t - tau));
                }
            }
            // Final point for trapezium rule - note first point is always 0
            if (t_idx >= 1)
            {
                I += AIF(t) / 2;
            }
            I = I * ktrans * m_dt;

            // VP component if vp>0
            f(t_idx + 1) = vp * AIF(t) + I;
        }
    }

    return f;
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
            Cp = m_ab * (1 - cos(m_mub * t)) + m_ab * m_ag * OrtonF(t, m_mug);
            c = vp * Cp;
            c += m_ab * m_ag * ktrans * (OrtonF(t, m_mug) + (((kep - m_mug) / m_ag) - 1) * OrtonF(t, kep)) / (kep - m_mug);
        }
        else
        {
            Cp = m_ab * m_ag * OrtonF(tb, m_mug) * exp(-m_mug * (t - tb));
            c = vp * Cp;
            c += m_ab * m_ag * ktrans * (OrtonF(tb, m_mug) * exp(-m_mug * (t - tb)) + (((kep - m_mug) / m_ag) - 1) * OrtonF(tb, kep) * exp(-kep * (t - tb))) / (kep - m_mug);
        }

        f(tp + 1) = c;
    }

    return f;
}

void DCEStdToftsFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Parameters that are inferred - extract and give sensible names
    int p = 1;
    double ktrans = params(p++);
    double kep, ve;
    if (m_infer_kep) {
        kep = params(p++);
        ve = ktrans / kep;
    }
    else {
        ve = params(p++);
        kep = ktrans / ve;
    }

    // Optional parameters to infer
    double vp = m_vp;
    if (m_infer_vp)
    {
        vp = params(p++);
    }

    // Standard DCE parameters - may be inferred
    double sig0 = m_sig0;
    double delay = m_delay;
    double t10 = m_t10;
    if (m_infer_sig0)
    {
        sig0 = params(p++);
    }
    if (m_infer_delay)
    {
        delay = params(p++);
    }
    if (m_infer_t10)
    {
        t10 = params(p++);
    }

    // Tissue concentration results
    ColumnVector C;
    if ((m_aif_type == "orton") && !m_force_conv)
    {
        C = GetConcentrationOrton(delay, vp, ktrans, kep);
    }
    else
    {
        C = GetConcentrationMeasuredAif(delay, vp, ktrans, kep);
    }

    // Convert concentration back to DCE signal
    result.ReSize(data.Nrows());
    for (int i = 1; i <= data.Nrows(); i++)
    {
        result(i) = SignalFromConcentration(C(i), t10, sig0);
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
