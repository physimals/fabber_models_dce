/*  fwdmodel_dce.h -Implements a convolution based model for DCE analysis

 Jesper Kallehauge, IBME

 Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */
#pragma once

#include "fabber_core/fwdmodel.h"

#include <newmat.h>

#include <string>
#include <vector>

/**
 * Implementation of the standard/extended Tofts model
 */
class DCEStdToftsFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    DCEStdToftsFwdModel()
    {
    }

    std::string GetDescription() const;
    virtual std::string ModelVersion() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;

    virtual void Initialize(FabberRunData &rundata);
    void GetParameterDefaults(std::vector<Parameter> &params) const;

    virtual void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;

private:
    // Mandatory
    double m_fa, m_tr, m_r1, m_dt;
    std::string m_aif_type;

    // Optional
    double m_vp, m_delay, m_t10, m_sig0;

    // Orton AIF
    double m_ab, m_ag, m_mub, m_mug;

    // Inference flags
    bool m_infer_vp, m_infer_delay, m_infer_t10, m_infer_sig0, m_infer_ve;

    // AIF as concentration curve
    NEWMAT::ColumnVector m_aif;

    double SignalFromConcentration(double C, double t10, double m0) const;
    double ConcentrationFromSignal(double s, double s0, double t10, double hct) const;
    NEWMAT::ColumnVector GetConcentrationMeasuredAif(double delay, double Vp, double Ktrans, double Kep) const;
    NEWMAT::ColumnVector aifshift(const NEWMAT::ColumnVector &aif, const double delay) const;
    NEWMAT::ColumnVector GetConcentrationOrton(double Vp, double Ktrans, double Ve) const;

    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DCEStdToftsFwdModel> registration;
};
