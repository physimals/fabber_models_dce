/**
 * fwdmodel_dce_tofts.h 
 *
 * Standard Tofts model for DCE analysis
 *
 * @author Martin Craig IBME Oxford
 *
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
 * Implementation of the standard/extended Tofts model
 */
class DCEStdToftsFwdModel : public DCEFwdModel
{
public:
    static FwdModel *NewInstance();

    DCEStdToftsFwdModel()
    {
    }

    std::string GetDescription() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;
    void Initialize(FabberRunData &rundata);
    void GetParameterDefaults(std::vector<Parameter> &params) const;
    void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;

private:
    // Initial values of model parameters - always inferred
    double m_ktrans, m_ve;

    // Optional model parameters - fixed when not inferred
    double m_vp;

    // Inference flags
    bool m_infer_vp, m_infer_kep;

    // Evaluation options
    bool m_force_conv;

    NEWMAT::ColumnVector GetConcentrationMeasuredAif(double delay, double Vp, double Ktrans, double Kep) const;
    NEWMAT::ColumnVector compute_tofts_model_measured_aif(double delay, double Vp, double Ktrans, double Ve) const;
    NEWMAT::ColumnVector GetConcentrationOrton(double delay, double Vp, double Ktrans, double Ve) const;

    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DCEStdToftsFwdModel> registration;
};
