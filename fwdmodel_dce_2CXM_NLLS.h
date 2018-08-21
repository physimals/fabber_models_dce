/*  fwdmodel_dce_2CXM_NLLS.h - Implementation of the non-linear least square solution of the two compartment exchange model

https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.25991

 Moss Zhao - IBME, Oxford

 Copyright (C) 2018 University of Oxford  */

/*  CCOPYRIGHT */
#pragma once

#include "fabber_core/fwdmodel.h"

#include <newmat.h>

#include <string>
#include <vector>

/**
 * Base class for DCE models as they share options
 */
class DCE2CXMNLLSFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    DCE2CXMNLLSFwdModel()
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
    double SignalFromConcentration(double C, double t10, double current_sig0) const;
    double ConcentrationFromSignal(double s, double s0, double t10, double hct) const;

private:
    // Mandatory
    double m_FA, m_TR, m_r1, m_dt;
    // AIF as concentration curve
    NEWMAT::ColumnVector m_aif;
    std::string m_aif_type, m_conv_method;

    // Optional
    double m_vp, m_delay, m_T10, m_sig0;

    // Orton AIF
    //double m_ab, m_ag, m_mub, m_mug;

    // Inference flags
    bool m_infer_vp, m_infer_delay, m_infer_t10, m_infer_sig0;

    // Other flags
    bool m_use_log;

    double LogOrNot(double p) const;
    NEWMAT::ColumnVector GetConcentrationMeasuredAif(double delay, double Vp, double Ktrans, double Kep) const;
    NEWMAT::ColumnVector aifshift(const NEWMAT::ColumnVector &aif, const double delay) const;
    NEWMAT::ColumnVector compute_concentration(const double Fp, const double PS, const double Vp, const double Ve, const NEWMAT::ColumnVector &aif) const;
    NEWMAT::ColumnVector compute_convolution_normal(const NEWMAT::ColumnVector &term_1, const NEWMAT::ColumnVector &term_2) const;
    NEWMAT::ColumnVector compute_convolution_iterative(const NEWMAT::ColumnVector &aif, const double T_term) const;
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DCE2CXMNLLSFwdModel> registration;
};
