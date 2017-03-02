/*  fwdmodel_asl_grase.h - Implements the GRASE model

 Michael Chappell, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce.h"

#include "fabber_core/fwdmodel.h"

#include <string>
#include <vector>

class DCE_2CXM_LLS_FwdModel : public DCEFwdModel
{
public:
    static FwdModel *NewInstance();

    // Virtual function overrides
    std::string GetDescription() const;
    virtual void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;
    virtual std::vector<std::string> GetUsage() const;
    virtual void NameParams(std::vector<std::string> &names) const;

    virtual void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DCE_2CXM_LLS_FwdModel> registration;
};
