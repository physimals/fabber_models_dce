/*  fwdmodel_dce_CTU_LLS.h - Implements the Linear Compartmental Uptake model

    Jesper Kallehauge, IBME

    Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabbercore/fwdmodel.h"
#include "fabbercore/inference.h"
#include <string>
using namespace std;

class DCE_CTU_LLS_FwdModel : public FwdModel {
public:
  static FwdModel* NewInstance();

  // Virtual function overrides
  virtual void Initialize(ArgsType& args);
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;
  virtual vector<string> GetUsage() const;
  virtual string ModelVersion() const;
                  
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const 
  { return 5 + (inferdelay?1:0); }

  virtual ~DCE_CTU_LLS_FwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

protected: 

  ColumnVector aifshift( const ColumnVector& aif, const float delta, const float hdelt ) const;
  void createconvmtx( LowerTriangularMatrix& A, const ColumnVector aifnew ) const;
  
// Constants

  // Lookup the starting indices of the parameters
  int Fp_index() const {return 1;}

  int Vp_index() const {  return 2;  }

  int PS_index() const {  return 3;  }

  int delta_index() const { return 3 + (inferdelay?1:0); }

  int T10_index() const { return 4 + (inferdelay?1:0); }

  int sig0_index() const { return 5 + (inferdelay?1:0); }

  //for ARD
  vector<int> ard_index;

    // scan parameters
  //double te;
  //double r2;
  double delt;
  double TR;
  double FA;
  double r1;
  double Tsat;
  double FA_radians;

  ColumnVector artsig;
  ColumnVector s;

  bool aifconc;

  bool inferdelay;
  bool doard;

  bool imageprior;

  string convmtx;
  string Acq_tech;
  
  private:
  /** Auto-register with forward model factory. */
  static FactoryRegistration<FwdModelFactory, DCE_CTU_LLS_FwdModel> registration;

};