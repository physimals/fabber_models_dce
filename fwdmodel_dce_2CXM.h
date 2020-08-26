/**
 * fwdmodel_dce_2CXM.h
 *
 * Implementation of the non-linear least square solution of the two compartment exchange model
 *
 * https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.25991
 *
 * @author Moss Zhao - IBME, Oxford
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
 * Two compartment exchange DCE forward model
 */
class DCE_2CXM_FwdModel : public DCEFwdModel
{
public:
    static FwdModel *NewInstance();

    DCE_2CXM_FwdModel()
    {
    }

    std::string GetDescription() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;
    void Initialize(FabberRunData &rundata);
    void GetParameterDefaults(std::vector<Parameter> &params) const;
    void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;

private:
    // Initial values of model parameters - always inferred
    double m_fp, m_ps, m_ve, m_vp;

    // Convolution method
    std::string m_conv_method;

    NEWMAT::ColumnVector compute_concentration(const double delay, const double Fp, const double PS, const double Vp, const double Ve) const;
    NEWMAT::ColumnVector compute_convolution_matrix(const double delay, const double T, const double T_plus, const double T_minus) const;
    NEWMAT::ColumnVector compute_convolution_iterative(const double delay, const double T_term) const;
    NEWMAT::ColumnVector compute_convolution_trap(const double delay, const double T, const double T_plus, const double T_minus) const;

    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DCE_2CXM_FwdModel> registration;
};
