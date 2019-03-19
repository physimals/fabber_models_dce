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

/**
 * Implementation of the standard/extended Tofts model
 */
class DCEFwdModel : public FwdModel
{
public:
    DCEFwdModel()
    {
    }

    virtual ~DCEFwdModel()
    {
    }

    virtual std::string ModelVersion() const;
    virtual void GetOptions(std::vector<OptionSpec> &opts) const;

    virtual void Initialize(FabberRunData &rundata);
    virtual void GetParameterDefaults(std::vector<Parameter> &params) const;
    
    virtual void InitVoxelPosterior(MVNDist &posterior) const;
    
protected:
    // Mandatory DCE configuration
    double m_fa, m_tr, m_r1, m_dt;
    std::string m_aif_type;

    // Standard DCE parameters (fixed when not inferred)
    double m_delay, m_t10, m_sig0;

    // Inference flags
    bool m_infer_delay, m_infer_t10, m_infer_sig0;

    // AIF as concentration curve
    NEWMAT::ColumnVector m_aif;

    // Other options
    mutable int m_sig0_idx;
    bool m_auto_init_delay;

    // Orton AIF parameters
    double m_ab, m_ag, m_mub, m_mug;

    double SignalFromConcentration(double C, double t10, double m0) const;
    double ConcentrationFromSignal(double s, double s0, double t10, double hct) const;
    NEWMAT::ColumnVector aifshift(const NEWMAT::ColumnVector &aif, const double delay) const;
    double OrtonF(double t, double a) const;
    double OrtonAIF(double t) const;
    double AIF(double t) const;
};

#if 0
class DCEFwdModel : public FwdModel
{
public:
    DCEFwdModel()
        : fp_idx(-1)
        , vp_idx(-1)
        , delta_idx(-1)
        , sig0_idx(-1)
        , t10_idx(-1)
        , ps_idx(-1)
        , ve_idx(-1)
        , ktrans_idx(-1)
    {
    }

    // Virtual function overrides
    void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string ModelVersion() const;
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

    virtual ~DCEFwdModel()
    {
        return;
    }

protected:
    void MakeParamIndex();
    std::map<std::string, int> paramIndex;

    NEWMAT::ColumnVector expConv(const NEWMAT::ColumnVector &aifnew, const float T,
        const NEWMAT::ColumnVector htsamp) const;
    NEWMAT::ColumnVector aifshift(const NEWMAT::ColumnVector &aif, const float delta,
        const float hdelt) const;
    void createconvmtx(NEWMAT::LowerTriangularMatrix &A, const NEWMAT::ColumnVector aifnew) const;

    // Indices of parameters
    int fp_idx, vp_idx, delta_idx, sig0_idx, t10_idx, ps_idx, ve_idx, ktrans_idx;

    double delt;
    double TR;
    double FA;
    double r1;
    double Tsat;
    double FA_radians;

    bool aifconc;

    bool inferdelay;
    bool doard;

    bool imageprior;

    std::string convmtx;
    std::string Acq_tech;

    NEWMAT::ColumnVector artsig;

    // for ARD
    std::vector<int> ard_index;
};

class DCEToftsFwdModel : public DCEFwdModel
{
public:
    static FwdModel *NewInstance();

    // Virtual function overrides
    std::string GetDescription() const;
    virtual void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;
    virtual std::vector<std::string> GetUsage() const;

    virtual ~DCEToftsFwdModel() { return; }
    virtual void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DCEToftsFwdModel> registration;
};
#endif
