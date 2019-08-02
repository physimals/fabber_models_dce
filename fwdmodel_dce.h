/**
 * fwdmodel_dce.h
 *
 * Base class for DCE models. Provides basic functions for dealing with AIF and
 * standard parameters (T1, delay, sig0).
 *
 * Martin Craig, IBME
 *
 * Copyright (C) 2016 University of Oxford  
 */

/*  CCOPYRIGHT */
#pragma once

#include "fabber_core/fwdmodel.h"

#include <newmat.h>

#include <string>
#include <vector>

/**
 * Base class for DCE models as they share options
 *
 * Should probably factor out AIF into a separate class one per implementation.
 * Currently have Orton (2008), Parker (2006) and user-specified (signal or concentration)
 */
class DCEFwdModel : public FwdModel
{
public:
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
    double OrtonF(double t, double a) const;
    double OrtonAIF(double t) const;
    double ParkerAIF(double t) const;
    double AIF(double t) const;
};
