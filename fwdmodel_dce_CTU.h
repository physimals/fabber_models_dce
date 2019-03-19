/**
 * fwdmodel_dce_ctu.h
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
#pragma once

#include "fwdmodel_dce.h"

#include <fabber_core/fwdmodel.h>

#include <newmat.h>

#include <string>
#include <vector>

/**
 * Compartment tissue update DCE model
 */
class DCE_CTU_FwdModel : public DCEFwdModel
{
public:
    static FwdModel *NewInstance();

    DCE_CTU_FwdModel()
    {
    }

    std::string GetDescription() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;
    void Initialize(FabberRunData &rundata);
    void GetParameterDefaults(std::vector<Parameter> &params) const;
    void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;
    
private:
    // Initial values of model parameters - always inferred
    double m_fp, m_ps, m_vp;

    // Convolution method
    std::string m_conv_method;

    NEWMAT::ColumnVector compute_concentration(double delay, double fp, double ps, double vp) const;
    NEWMAT::ColumnVector compute_convolution_matrix(double delay, const NEWMAT::ColumnVector &term_2) const;
    NEWMAT::ColumnVector compute_convolution_iterative(const double delay, const double fp, const double ktrans, const double Tp) const;
    NEWMAT::ColumnVector compute_convolution_trap(const double delay, const double fp, const double ktrans, const double Tp) const;

    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DCE_CTU_FwdModel> registration;
};
