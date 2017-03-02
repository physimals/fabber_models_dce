/*  fwdmodel_dce_CTU_LLS.h - Implements the Linear Compartmental Uptake model

 Jesper Kallehauge, IBME

 Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce.h"

#include "fabber_core/fwdmodel.h"

#include <string>
#include <vector>

class DCE_CTU_LLS_FwdModel : public DCEFwdModel
{
public:
    static FwdModel *NewInstance();

    // Virtual function overrides
    std::string GetDescription() const;
    virtual void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;
    virtual std::vector<std::string> GetUsage() const;

    virtual void NameParams(std::vector<std::string> &names) const;
    virtual void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const;

protected:
private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DCE_CTU_LLS_FwdModel> registration;
};
