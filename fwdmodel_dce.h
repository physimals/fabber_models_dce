/*  fwdmodel_dce.h -Implements a convolution based model for DCE analysis

 Jesper Kallehauge, IBME

 Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */
#pragma once

#include "fabber_core/fwdmodel.h"
#include "fabber_core/inference.h"
#include <string>
using namespace std;

using namespace NEWMAT;

/**
 * Base class for DCE models as they share options
 */
class DCEFwdModel: public FwdModel
{
public:
	DCEFwdModel() :
			fp_idx(-1), vp_idx(-1), delta_idx(-1), sig0_idx(-1), t10_idx(-1), ps_idx(-1), ve_idx(-1), ktrans_idx(-1)
	{
	}
	// Virtual function overrides
	void GetOptions(vector<OptionSpec> &opts) const;
	virtual string ModelVersion() const;
	virtual void Initialize(ArgsType& args);

	virtual void DumpParameters(const ColumnVector& vec, const string& indents = "") const;

	virtual void NameParams(vector<string>& names) const;
	virtual int NumParams() const
	{
		vector<string> names;
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

	ColumnVector expConv(const ColumnVector& aifnew, const float T, const ColumnVector htsamp) const;
	ColumnVector aifshift(const ColumnVector& aif, const float delta, const float hdelt) const;
	void createconvmtx(LowerTriangularMatrix& A, const ColumnVector aifnew) const;

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

	string convmtx;
	string Acq_tech;

	ColumnVector artsig;

	//for ARD
	vector<int> ard_index;
};

class DCEToftsFwdModel: public DCEFwdModel
{
public:
	static FwdModel* NewInstance();

	// Virtual function overrides
	std::string GetDescription() const;
	virtual void Evaluate(const ColumnVector& params, ColumnVector& result) const;
	virtual vector<string> GetUsage() const;

	virtual ~DCEToftsFwdModel()
	{
		return;
	}

	virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

private:
	/** Auto-register with forward model factory. */
	static FactoryRegistration<FwdModelFactory, DCEToftsFwdModel> registration;

};
