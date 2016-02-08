/*  fwdmodel_dce_ETM.h - Implements the Extended Tofts model

    Jesper Kallehauge, IBME

    Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabbercore/fwdmodel.h"
#include "fabbercore/inference.h"
#include <string>
using namespace std;

class DCE_ETM_FwdModel : public FwdModel {
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

  virtual ~DCE_ETM_FwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

protected: 

  ColumnVector aifshift( const ColumnVector& aif, const float delta, const float hdelt ) const;
  ColumnVector expConv( const ColumnVector& aifnew,  const float T, const ColumnVector htsamp) const;
  void createconvmtx( LowerTriangularMatrix& A, const ColumnVector aifnew ) const;
  
// Constants

  // Lookup the starting indices of the parameters
  int Ktrans_index() const {return 1;}

  int Ve_index() const {  return 2;  }

  int Vp_index() const {  return 3;  }

  int delta_index() const { return 3 + (inferdelay?1:0); }

  int sig0_index() const { return 4 + (inferdelay?1:0); }
  
  int T10_index() const { return 5 + (inferdelay?1:0); }


  //for ARD
  vector<int> ard_index;

  // scan parameters
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

  string Acq_tech;
  string convmtx;
  
  private:
  /** Auto-register with forward model factory. */
  static FactoryRegistration<FwdModelFactory, DCE_ETM_FwdModel> registration;

};