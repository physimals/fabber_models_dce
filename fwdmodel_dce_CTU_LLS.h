/*  fwdmodel_dce_CTU_LLS.h - Implements the Linear Compartmental Uptake model

 Jesper Kallehauge, IBME

 Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce.h"

#include "fabber_core/fwdmodel.h"
#include "fabber_core/inference.h"
#include <string>
using namespace std;

using namespace NEWMAT;

class DCE_CTU_LLS_FwdModel: public DCEFwdModel
{
public:
	static FwdModel* NewInstance();

	// Virtual function overrides
	std::string GetDescription() const;
	virtual void Evaluate(const ColumnVector& params, ColumnVector& result) const;
	virtual vector<string> GetUsage() const;

	virtual void NameParams(vector<string>& names) const;

	virtual ~DCE_CTU_LLS_FwdModel()
	{
		return;
	}

	virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

protected:

private:
	/** Auto-register with forward model factory. */
	static FactoryRegistration<FwdModelFactory, DCE_CTU_LLS_FwdModel> registration;

};
