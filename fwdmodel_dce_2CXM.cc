/**
 * fwdmodel_dce_2CXM.cc
 *
 * Implementation of the non-linear least square solution of the two compartment exchange model
 *
 * https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.25991
 * 
 * Moss Zhao - IBME, Oxford
 *
 * Copyright (C) 2018 University of Oxford
 */

/*  CCOPYRIGHT */

#include "fwdmodel_dce_2CXM.h"

#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>

#include <miscmaths/miscprob.h>
#include <newimage/newimageall.h>

#include <newmatio.h>

#include <iostream>
#include <stdexcept>

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DCE_2CXM_FwdModel> DCE_2CXM_FwdModel::registration("dce_2CXM");

std::string DCE_2CXM_FwdModel::GetDescription() const
{
    return "Non-linear least square solution of the two compartment exchange model";
}

static OptionSpec OPTIONS[] = {
    { "fp", OPT_FLOAT, "Flow in min-1", OPT_NONREQ, "0.3" },
    { "ps", OPT_FLOAT, "Permeability surface area product in min-1", OPT_NONREQ, "0.3" },
    { "vp", OPT_FLOAT, "Plasma volume in decimal between zero and one", OPT_NONREQ, "0.3" },
    { "ve", OPT_FLOAT, "Extracellular space volume in decimal between zero and one", OPT_NONREQ, "0.3" },
    { "convolution-method", OPT_STR, "Method to compute convolution, normal or iterative. Default is iterative", OPT_REQ, "" },
    { "" },
};

void DCE_2CXM_FwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    DCEFwdModel::GetOptions(opts);
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

void DCE_2CXM_FwdModel::Initialize(FabberRunData &rundata)
{
    DCEFwdModel::Initialize(rundata);
    
    // Initial values of main parameters
    m_fp = rundata.GetDoubleDefault("fp", 0.3);
    m_ps = rundata.GetDoubleDefault("ps", 0.3);
    m_vp = rundata.GetDoubleDefault("vp", 0.3);
    m_ve = rundata.GetDoubleDefault("ve", 0.3);

    // Other model options
    m_conv_method = rundata.GetStringDefault("convolution-method", "normal");
}

void DCE_2CXM_FwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    // Basic model parameters
    int p=0;
    params.push_back(Parameter(p++, "fp", DistParams(m_fp, 100), DistParams(m_fp, 100), PRIOR_NORMAL, TRANSFORM_ABS()));
    params.push_back(Parameter(p++, "ps", DistParams(m_ps, 100), DistParams(m_ps, 100), PRIOR_NORMAL, TRANSFORM_ABS()));
    params.push_back(Parameter(p++, "ve", DistParams(m_ve, 10), DistParams(m_ve, 10), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    params.push_back(Parameter(p++, "vp", DistParams(m_vp, 10), DistParams(m_vp, 10), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    
    // Standard DCE parameters
    DCEFwdModel::GetParameterDefaults(params);
}

ColumnVector DCE_2CXM_FwdModel::compute_concentration(const double delay, const double fp, const double ps, const double vp, const double ve) const
{
    // Equation 1 of the paper
    double Tp = vp / fp;
    double Te = ve / ps;
    double T = (vp + ve) / fp;

    // Equation 9 of the paper
    double T_plus = 0.5 * (T + Te + sqrt(pow((T + Te), 2.0) - 4 * Tp * Te));
    double T_minus = 0.5 * (T + Te - sqrt(pow((T + Te), 2.0) - 4 * Tp * Te));

    // Result concentration
    ColumnVector current_concentration;

    if(m_conv_method == "normal") {
        // Compute convolution using normal technique
        
        ColumnVector convolution_result = compute_convolution_normal(delay, T, T_plus, T_minus);
        //cout << convolution_result.t() << endl;
        current_concentration = fp * convolution_result;
    }
    else if(m_conv_method == "msc") {
        // Compute convolution using msc technique
        
        ColumnVector convolution_result = compute_convolution_msc(delay, T, T_plus, T_minus);
        cout << convolution_result.t() << endl;
        current_concentration = fp * convolution_result;
    }
    else if(m_conv_method == "iterative") {
        // Compute convolution using iterative technique 
        // We use the method in the appendix of this paper to compute the convolution
        //Some properties of convolution. X represents convolution operation
        //f X g = g X f
        //a * (f X g) = (a * f X g)
        //(f + h) X g = f X g + h X g
        //We need to rearrange the terms in Equation 7 using these properties
        ColumnVector sum_1 = T_plus * (T - T_minus) * compute_convolution_iterative(delay, T_plus) / (T_plus - T_minus);
        ColumnVector sum_2 = T_minus * (T_plus - T) * compute_convolution_iterative(delay, T_minus) / (T_plus - T_minus);

        // Compute concentration
        current_concentration = fp * (sum_1 + sum_2);
    }

    return current_concentration;
}

// Compute convolution using normal method
// Assuming that the two input vectors have the same length
// https://stackoverflow.com/questions/24518989/how-to-perform-1-dimensional-valid-convolution
ColumnVector DCE_2CXM_FwdModel::compute_convolution_normal(const double delay, const double T, const double T_plus, const double T_minus) const
{
    // Get possibly shifted AIF
    ColumnVector aif(data.Nrows());
    for (int t_idx=0; t_idx<data.Nrows(); t_idx++) {
        double t = double(t_idx) * m_dt - delay;
        aif(t_idx+1) = AIF(t);
    }

    // We only need the first Nrows() terms
    ColumnVector convolution_result(data.Nrows());

    // Bracket term in Equation 7
    ColumnVector bracket_term(data.Nrows());
    for (int t_index = 0; t_index < data.Nrows(); t_index++) {
        double t = t_index * m_dt - delay;
        bracket_term(t_index + 1) = (T - T_minus) * (exp((-1) * t / T_plus)) / (T_plus - T_minus) + (T_plus - T) * (exp((-1) * t / T_minus)) / (T_plus - T_minus);
    }

    for (int t_index = 0; t_index < aif.Nrows(); t_index++) {
        convolution_result(t_index + 1) = 0;
        int jmn = (t_index >= aif.Nrows() - 1)? t_index - (aif.Nrows() - 1) : 0;
        int jmx = (t_index <  aif.Nrows() - 1)? t_index                             : aif.Nrows() - 1;
        for (int j = jmn; j <= jmx; j++) {
            convolution_result(t_index + 1) += (aif(j + 1) * bracket_term(t_index - j + 1));
        }
    }
    return convolution_result;
}

// Compute convolution using trapezium rule
ColumnVector DCE_2CXM_FwdModel::compute_convolution_msc(const double delay, const double T, const double T_plus, const double T_minus) const
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
                    // Bracket term in Equation 7 from paper
                    double bracket_term = (T - T_minus) * (exp(-(t - tau) / T_plus)) / (T_plus - T_minus) + (T_plus - T) * (exp(-(t - tau) / T_minus)) / (T_plus - T_minus);
                    I += AIF(tau) * bracket_term;
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

// Implementation of appendix of the paper
ColumnVector DCE_2CXM_FwdModel::compute_convolution_iterative(const double delay, const double T_term) const
{
    // Get possibly shifted AIF
    ColumnVector aif(data.Nrows());
    for (int t_idx=0; t_idx<data.Nrows(); t_idx++) {
        double t = double(t_idx) * m_dt - delay;
        aif(t_idx+1) = AIF(t);
    }
    
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

void DCE_2CXM_FwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Parameters that are inferred - extract and give sensible names
    int p = 1;
    double fp = params(p++);
    double ps = params(p++);
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

    // Tissue concentration results
    ColumnVector concentration_tissue = compute_concentration(delay, fp, ps, vp, ve);
    
    // Converts concentration back to DCE signal
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
            break;
        }
    }
}

FwdModel *DCE_2CXM_FwdModel::NewInstance()
{
    return new DCE_2CXM_FwdModel();
}
