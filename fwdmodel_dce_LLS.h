/*  fwdmodel_asl_grase.h - Implements the GRASE model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce.h"

#include "fabber_core/fwdmodel.h"
#include "fabber_core/inference.h"
#include <string>
using namespace std;

using namespace NEWMAT;

class DCE_LLS_FwdModel : public DCEFwdModel {
public:
  static FwdModel* NewInstance();

  // Virtual function overrides
  std::string GetDescription() const;
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;
  virtual vector<string> GetUsage() const;
  virtual string ModelVersion() const;
                                               
  virtual void NameParams(vector<string>& names) const;     
  
  virtual ~DCE_LLS_FwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

protected: 

  private:
  /** Auto-register with forward model factory. */
  static FactoryRegistration<FwdModelFactory, DCE_LLS_FwdModel> registration;

};
