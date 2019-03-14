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
 * Base class for DCE models as they share options
 */
class DCEStdToftsFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    DCEStdToftsFwdModel()
    {
    }

    void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string ModelVersion() const;
    std::string GetDescription() const;
    virtual std::vector<std::string> GetUsage() const;

    virtual void Initialize(ArgsType &args);

    virtual void DumpParameters(const NEWMAT::ColumnVector &vec,
        const std::string &indents = "") const;

    virtual void NameParams(std::vector<std::string> &names) const;
    virtual int NumParams() const
    {
        std::vector<std::string> names;
        NameParams(names);
        return names.size();
    }

    virtual void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const;

    //virtual bool Gradient(const NEWMAT::ColumnVector &params, NEWMAT::Matrix &grad) const;
    virtual void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;
    double SignalFromConcentration(double C, double t10, double m0) const;
    double ConcentrationFromSignal(double s, double s0, double t10, double hct) const;

private:
    // Mandatory
    double m_FA, m_TR, m_r1, m_dt;
    std::string m_aif_type;
    // Optional initial values
    double initial_Ktrans, initial_Vp, initial_Ve;

    // Optional
    double m_vp, m_delay, m_T10, m_sig0;

    // Orton AIF
    double m_ab, m_ag, m_mub, m_mug;

    // Inference flags
    bool m_infer_vp, m_infer_delay, m_infer_t10, m_infer_sig0;

    // Other flags
    bool m_use_log;

    // AIF as concentration curve
    NEWMAT::ColumnVector m_aif;

    double LogOrNot(double p) const;
    NEWMAT::ColumnVector GetConcentrationMeasuredAif(double delay, double Vp, double Ktrans, double Kep) const;
    NEWMAT::ColumnVector compute_tofts_model_measured_aif(double delay, double Vp, double Ktrans, double Ve) const;
    NEWMAT::ColumnVector aifshift(const NEWMAT::ColumnVector &aif, const double delay) const;
    NEWMAT::ColumnVector GetConcentrationOrton(double Vp, double Ktrans, double Ve) const;

    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DCEStdToftsFwdModel> registration;
};
