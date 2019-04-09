/**
 * fwdmodel_dce_ctu.cc
 *
 * Implementation of the compartment tissue update model
 *
 * https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.26324
 *
 * Moss Zhao - IBME, Oxford
 *
 * Copyright (C) 2018 University of Oxford  
 */

/*  CCOPYRIGHT */

#include "fwdmodel_dce_CTU.h"

#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>

#include <miscmaths/miscprob.h>
#include <newimage/newimageall.h>

#include <newmatio.h>

#include <iostream>
#include <stdexcept>
#include <math.h>

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DCE_CTU_FwdModel> DCE_CTU_FwdModel::registration("dce_CTU");

std::string DCE_CTU_FwdModel::GetDescription() const
{
    return "Compartment tissue update model";
}

static OptionSpec OPTIONS[] = {
    { "fp", OPT_FLOAT, "Prior value for flow in min-1", OPT_NONREQ, "0.3" },
    { "ps", OPT_FLOAT, "Prior value for permeability surface area product in min-1", OPT_NONREQ, "0.3" },
    { "vp", OPT_FLOAT, "Prior value for plasma volume in decimal between zero and one", OPT_NONREQ, "0.3" },
    { "conv-method", OPT_STR, "Method to compute convolution, trapezium, matrix or iterative.", OPT_REQ, "trapezium" },
    { "" },
};

void DCE_CTU_FwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    DCEFwdModel::GetOptions(opts);
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

void DCE_CTU_FwdModel::Initialize(FabberRunData &rundata)
{
    DCEFwdModel::Initialize(rundata);
    
    // Initial values of the parameters
    m_fp = rundata.GetDoubleDefault("fp", 0.5);
    m_ps = rundata.GetDoubleDefault("ps", 0.5);
    m_vp = rundata.GetDoubleDefault("vp", 0.05);

    // Other model options
    m_conv_method = rundata.GetStringDefault("conv-method", "trapezium");
}

void DCE_CTU_FwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    // Basic model parameters
    int p=0;
    params.push_back(Parameter(p++, "fp", DistParams(m_fp, 1e5), DistParams(m_fp, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
    params.push_back(Parameter(p++, "ps", DistParams(m_ps, 1e5), DistParams(m_ps, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
    params.push_back(Parameter(p++, "vp", DistParams(m_vp, 1), DistParams(m_vp, 1), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    
    // Standard DCE parameters
    DCEFwdModel::GetParameterDefaults(params);
}

// Implementation of Equation 4 of the paper
ColumnVector DCE_CTU_FwdModel::compute_concentration(double delay, double fp, double ps, double vp) const
{
    if(fp == 0) {
        fp = 0.00001;
    }
    if(ps == 0) {
        ps = 0.00001;
    }
    if(vp == 0) {
        vp = 0.00001;
    }

    double Tp = vp / (fp + ps);
    double E = ps / (ps + fp);
    double ktrans = E * fp;

    // Result concentration
    ColumnVector current_concentration(data.Nrows());

    if(m_conv_method == "matrix") {
        // Compute convolution using matrix multiplication
        ColumnVector convolution_result(data.Nrows());

        // Bracket term in Equation 4
        ColumnVector bracket_term(data.Nrows());

        for (int t_index = 0; t_index < data.Nrows(); t_index++) {
            double current_t_value = t_index * m_dt - delay;
            bracket_term(t_index + 1) = fp * exp(-current_t_value / Tp) + ktrans * (1 - exp(-current_t_value / Tp));
            current_concentration(t_index + 1) = 0.0; // Initialize the concentration result vector. It seems that C++ likes this. :O
        }

        // Do convolution. 
        current_concentration = m_dt * compute_convolution_matrix(delay, bracket_term);
    }
    else if(m_conv_method == "trapezium") {
        // Compute convolution using trapezium rule
        current_concentration = compute_convolution_trap(delay, fp, ktrans, Tp);
    }
    else if(m_conv_method == "iterative") {
        // Need to work on this (you need to re-write the equation)
        // Compute convolution using iterative technique 
        // We use the method in the appendix of this paper to compute the convolution
        //Some properties of convolution. X represents convolution operation
        //f X g = g X f
        //a * (f X g) = (a * f X g)
        //(f + h) X g = f X g + h X g
        //We need to rearrange the terms in Equation 4 using these properties
        //"""We need to rearrange Equation 4 of the model paper into the following way:"""
        //"""C(t) = (fp - ktrans) * Tp * (exp(-t/Tp))/Tp X Ca + ktrans * Ca"""
        current_concentration = compute_convolution_iterative(delay, fp, ktrans, Tp);
    }

    return current_concentration;
}

// Compute convolution using trapezium rule
ColumnVector DCE_CTU_FwdModel::compute_convolution_trap(const double delay, const double fp, const double ktrans, const double Tp) const
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
                    // Bracket term in Equation 4 from paper
                    double bracket_term = fp * exp(-(t-tau) / Tp) + ktrans * (1 - exp(-(t-tau) / Tp));
                    I += AIF(tau) * bracket_term;
                }
            }
            if (t_idx >= 1) 
            {
                // Last point half contribution in trapezium rule. note first point
                // is always zero since AIF(0) = 0
                I += fp * AIF(t) / 2;
            }
            convolution_result(t_idx + 1) = I * m_dt;
        }
    }

    return convolution_result;
}

// Compute convolution using matrix multiplication
// Assuming that the two input vectors have the same length
// https://stackoverflow.com/questions/24518989/how-to-perform-1-dimensional-valid-convolution
ColumnVector DCE_CTU_FwdModel::compute_convolution_matrix(double delay, const NEWMAT::ColumnVector &term_2) const
{
    // Get possibly shifted AIF
    ColumnVector aif(data.Nrows());
    for (int t_idx=0; t_idx<data.Nrows(); t_idx++) {
        double t = double(t_idx) * m_dt - delay;
        aif(t_idx+1) = AIF(t);
    }

    // We only need the first Nrows() terms
    ColumnVector convolution_result(aif.Nrows());

    for (int t_index = 0; t_index < aif.Nrows(); t_index++) {
        convolution_result(t_index + 1) = 0;
        int jmn = (t_index >= aif.Nrows() - 1)? t_index - (aif.Nrows() - 1) : 0;
        int jmx = (t_index <  aif.Nrows() - 1)? t_index                      : aif.Nrows() - 1;
        for (int j = jmn; j <= jmx; j++) {
            convolution_result(t_index + 1) += (aif(j + 1) * term_2(t_index - j + 1));
        }
    }
    return convolution_result;
}

// Implementation of appendix of the paper
ColumnVector DCE_CTU_FwdModel::compute_convolution_iterative(const double delay, const double fp, const double ktrans, const double Tp) const
{
    // Get possibly shifted AIF
    ColumnVector aif(data.Nrows());
    for (int t_idx=0; t_idx<data.Nrows(); t_idx++) {
        double t = double(t_idx) * m_dt - delay;
        aif(t_idx+1) = AIF(t);
    }
    
    int nhtpts = aif.Nrows();
    ColumnVector f_vector(nhtpts);

    if(Tp == 0) {
        f_vector = aif;
    }
    else {
        ColumnVector E0_vector(nhtpts - 1);
        ColumnVector E1_vector(nhtpts - 1);
        ColumnVector xi_vector(nhtpts - 1);
        ColumnVector aif_prime_vector(nhtpts - 1);

        for (int i = 2; i <= nhtpts; i++) {
            xi_vector(i - 1) = m_dt / Tp; // Because in the equation (A5) t_i+1 - t_i is always m_dt
            aif_prime_vector(i - 1) = (aif(i) - aif(i - 1)) / m_dt;
            E0_vector(i - 1) = 1 - exp(-xi_vector(i - 1));
            E1_vector(i - 1) = xi_vector(i - 1) - E0_vector(i - 1);
        }

        f_vector(1) = 0.0;

        for (int i = 2; i <= nhtpts; i++) {
            f_vector(i) = exp(-xi_vector(i - 1)) * f_vector(i - 1) + aif(i - 1) * E0_vector(i - 1) + aif_prime_vector(i - 1) * Tp * E1_vector(i - 1);
        }

    }

    return (fp - ktrans) * Tp * f_vector + ktrans * aif;
}

void DCE_CTU_FwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Parameters that are inferred - extract and give sensible names
    int p = 1;
    double fp = params(p++);
    double ps = params(p++);
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

    ColumnVector concentration_tissue = compute_concentration(delay, fp, ps, vp);

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

FwdModel *DCE_CTU_FwdModel::NewInstance()
{
    return new DCE_CTU_FwdModel();
}
