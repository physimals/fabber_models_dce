/**
 * fwdmodel_dce_AATH.h
 *
 * An Adiabatic Approximation to the Tissue Homogeneity Model for Water Exchange in the 
 * Brain: I. Theoretical Derivation
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
 * Adiabatic Approximation to the Tissue Homogeneity DCE forward model
 */
class DCE_AATH_FwdModel : public DCEFwdModel
{
public:
    static FwdModel *NewInstance();

    DCE_AATH_FwdModel()
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

    NEWMAT::ColumnVector compute_concentration(const double delay, const double Fp, const double PS, const double Vp, const double Ve) const;
    NEWMAT::ColumnVector compute_convolution_matrix(const NEWMAT::ColumnVector &term_1, const NEWMAT::ColumnVector &term_2) const;
    NEWMAT::ColumnVector compute_convolution_trap(const double delay, const double k_adb) const;

    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DCE_AATH_FwdModel> registration;
};
