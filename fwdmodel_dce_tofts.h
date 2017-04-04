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
    virtual void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;
    float SignalFromConcentration(float C, float t10,float m0) const;
    float ConcentrationFromSignal(float s, float s0, float t10, float hct) const;
private:
    double m_FA, m_TR, m_r1, m_dt, m_vp, m_delay, m_T10, m_sig0;
    bool m_infer_vp, m_infer_delay, m_infer_t10, m_infer_sig0, m_use_log;
    NEWMAT::ColumnVector aif;
    
    
    float LogOrNot(float p) const;
    NEWMAT::ColumnVector convolution(const NEWMAT::ColumnVector &myaif, float Vp, float Ktrans, float Kep) const;
    NEWMAT::ColumnVector aifshift(const NEWMAT::ColumnVector &aif, const float delay) const;
    
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DCEStdToftsFwdModel> registration;
};
