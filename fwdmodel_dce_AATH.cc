/**
 * fwdmodel_dce_AATH.cc
 *
 * An Adiabatic Approximation to the Tissue Homogeneity Model for Water Exchange in the 
 * Brain: I. Theoretical Derivation
 *
 * https://journals.sagepub.com/doi/10.1097/00004647-199812000-00011
 *
 * Moss Zhao - IBME, Oxford
 *
 * Copyright (C) 2018 University of Oxford  
 */

/*  CCOPYRIGHT */

#include "fwdmodel_dce_AATH.h"

#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>

#include <miscmaths/miscprob.h>
#include <newimage/newimageall.h>

#include <newmatio.h>

#include <iostream>
#include <stdexcept>

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DCE_AATH_FwdModel> DCE_AATH_FwdModel::registration("dce_AATH");

std::string DCE_AATH_FwdModel::GetDescription() const
{
    return "Adiabatic Approximation to the Tissue Homogeneity Model";
}

static OptionSpec OPTIONS[] = {
    { "infer-fp_AATH", OPT_BOOL, "Infer Fp in AATH model.", OPT_NONREQ, "" },
    { "fp", OPT_FLOAT, "Flow in min-1", OPT_NONREQ, "0.3" },
    { "infer-ps_AATH", OPT_BOOL, "Infer PS in AATH model.", OPT_NONREQ, "" },
    { "ps", OPT_FLOAT, "Permeability surface area product in min-1", OPT_NONREQ, "0.3" },
    { "vp", OPT_FLOAT, "Plasma volume in decimal between zero and one", OPT_NONREQ, "0.3" },
    { "ve", OPT_FLOAT, "Extracellular space volume in decimal between zero and one", OPT_NONREQ, "0.3" },
    { "" },
};

void DCE_AATH_FwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

void DCE_AATH_FwdModel::Initialize(FabberRunData &rundata)
{
    DCEFwdModel::Initialize(rundata);
    
    // Initial values of main parameters
    m_fp = rundata.GetDoubleDefault("fp", 0.3);
    m_ps = rundata.GetDoubleDefault("ps", 0.3);
    m_vp = rundata.GetDoubleDefault("vp", 0.3);
    m_ve = rundata.GetDoubleDefault("ve", 0.3);

    // ps and fp are AATH model parameters. At least one of them can be specified.
    // Both of them can be estimated at the same time but will not be accurate due to a simple multiplication problem
    // a * b = c. We can't estimate a or b by just knowing c.
    m_infer_fp_AATH = rundata.ReadBool("infer-fp_AATH");
    m_infer_ps_AATH = rundata.ReadBool("infer-ps_AATH");
}

void DCE_AATH_FwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    // Basic model parameters
    int p=0;
    if (m_infer_fp_AATH) {
        params.push_back(Parameter(p++, "fp", DistParams(m_fp, 100), DistParams(m_fp, 100), PRIOR_NORMAL, TRANSFORM_ABS()));
    }
    if (m_infer_ps_AATH) {
        params.push_back(Parameter(p++, "ps", DistParams(m_ps, 100), DistParams(m_ps, 100), PRIOR_NORMAL, TRANSFORM_ABS()));
    }
    params.push_back(Parameter(p++, "ve", DistParams(m_ve, 10), DistParams(m_ve, 10), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    params.push_back(Parameter(p++, "vp", DistParams(m_vp, 10), DistParams(m_vp, 10), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    
    // Standard DCE parameters
    DCEFwdModel::GetParameterDefaults(params);
}

ColumnVector DCE_AATH_FwdModel::compute_concentration(double delay, double fp, double ps, double vp, double ve) const
{   
    if(ve == 0) {
        ve = 0.0001;
    }
    if(ve == 1) {
        ve = 0.9999;
    }
    if(fp == 0) {
        fp = 0.0001;
    }

    double E = 1 - exp(-ps / fp);
    double k_adb = E * fp / ve;
    
    // Get possibly shifted AIF
    ColumnVector aif(data.Nrows());
    for (int t_idx=0; t_idx<data.Nrows(); t_idx++) {
        double t = double(t_idx) * m_dt - delay;
        aif(t_idx+1) = AIF(t);
    }

    // Do convolution.
    //ColumnVector exp_term(data.Nrows());
    //for (int t_index = 0; t_index < data.Nrows(); t_index++) {
    //    double current_t_value = t_index * m_dt - m_delay;
    //    exp_term(t_index + 1) = exp(-k_adb * current_t_value);
    //}
    //ColumnVector convolution_result = compute_convolution_matrix(aif, exp_term);
    ColumnVector convolution_result = compute_convolution_trap(delay, k_adb);
    
    return vp * aif + E * fp * convolution_result;
}

// Compute convolution using matrix multiplication
// Assuming that the two input vectors have the same length
// https://stackoverflow.com/questions/24518989/how-to-perform-1-dimensional-valid-convolution
ColumnVector DCE_AATH_FwdModel::compute_convolution_matrix(const NEWMAT::ColumnVector &term_1, const NEWMAT::ColumnVector &term_2) const
{
    ColumnVector convolution_result(term_1.Nrows());

    for (int t_index = 0; t_index < term_1.Nrows(); t_index++) {
        convolution_result(t_index + 1) = 0;
        int jmn = (t_index >= term_1.Nrows() - 1)? t_index - (term_1.Nrows() - 1) : 0;
        int jmx = (t_index <  term_1.Nrows() - 1)? t_index                      : term_1.Nrows() - 1;
        for (int j = jmn; j <= jmx; j++) {
            convolution_result(t_index + 1) += (term_1(j + 1) * term_2(t_index - j + 1));
        }
    }
    return m_dt * convolution_result;
}

// Compute convolution using trapezium rule
ColumnVector DCE_AATH_FwdModel::compute_convolution_trap(const double delay, const double k_adb) const
{
    ColumnVector convolution_result(data.Nrows());

    for (int t_idx = 0; t_idx < data.Nrows(); t_idx++)
    {
        double t = t_idx * m_dt - delay;
        if (t <= 0) 
        {
            // AIF strictly zero at t <= 0
            convolution_result(t_idx + 1) = 0;
        }
        else 
        {
            // Convolution integral: INTEGRAL 0->t {AIF(t) * bracket_term(t-tau)} d tau
            double I = 0;
            for (int tau_idx = 0; tau_idx < t_idx; tau_idx++)
            {
                double tau = tau_idx * m_dt - delay;
                if (tau > 0) {
                    I += AIF(tau) * exp(-k_adb * (t-tau));
                }
            }
            if (t_idx >= 1) 
            {
                // Last point half contribution in trapezium rule. note first point
                // is always zero since AIF(0) = 0
                I += AIF(t) / 2;
            }
            convolution_result(t_idx + 1) = I * m_dt;
        }
    }

    return convolution_result;
}

void DCE_AATH_FwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Parameters that are inferred - extract and give sensible names
    int p = 1;
    double fp = m_fp;
    double ps = m_ps;
    if (m_infer_fp_AATH) {
        fp = params(p++);
    }
    if (m_infer_ps_AATH) {
        ps = params(p++);
    }
    double ve = params(p++);
    double vp = params(p++);

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

    ColumnVector concentration_tissue = compute_concentration(delay, fp, ps, vp, ve);

    // Converts concentration back to MRI signal
    result.ReSize(data.Nrows());
    for (int i = 1; i <= data.Nrows(); i++)
    {   
        result(i) = SignalFromConcentration(concentration_tissue(i), t10, sig0);
    }

    for (int i = 1; i <= data.Nrows(); i++)
    {
        if (isnan(result(i)) || isinf(result(i)))
        {
            LOG << "Warning NaN or inf in result" << endl;
            LOG << "result: " << result.t() << endl;
            LOG << "params: " << params.t() << endl;

            result = 0.0;
            //getchar();
            break;

        }
    }
}

FwdModel *DCE_AATH_FwdModel::NewInstance()
{
    return new DCE_AATH_FwdModel();
}
