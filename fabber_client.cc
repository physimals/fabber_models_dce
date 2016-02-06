/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

    Michael Chappell, FMRIB Image Analysis & IBME QuBIc Groups

    Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabbercore/fabber_core.h"

// DCE models to be included from library
#include "fwdmodel_dce.h"
#include "fwdmodel_dce_LLS.h"
#include "fwdmodel_dce_Patlak.h"
#include "fwdmodel_dce_ETM.h"
#include "fwdmodel_dce_ETM_LLS.h"
#include "fwdmodel_dce_CTU.h"
#include "fwdmodel_dce_CTU_LLS.h"
#include "fwdmodel_dce_2CXM.h"
#include "fwdmodel_dce_2CXM_LLS.h"
#include "fwdmodel_dce_AATH.h"



int main(int argc, char** argv) {

  //add the DCE models - these will autoregister at this point
  DCEFwdModel::NewInstance();
  DCE_LLS_FwdModel::NewInstance();
  DCE_Patlak_FwdModel::NewInstance();
  DCE_ETM_FwdModel::NewInstance();
  DCE_ETM_LLS_FwdModel::NewInstance();
  DCE_CTU_FwdModel::NewInstance();
  DCE_CTU_LLS_FwdModel::NewInstance();
  DCE_2CXM_FwdModel::NewInstance();
  DCE_2CXM_LLS_FwdModel::NewInstance();
  DCE_AATH_FwdModel::NewInstance();

  
  return execute(argc, argv);

}