/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

    Michael Chappell, FMRIB Image Analysis & IBME QuBIc Groups

    Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabbercore/fabber_core.h"

// CEST models to be included from library
#include "fwdmodel_dce.h"
#include "fwdmodel_dce_ETM.h"
#include "fwdmodel_dce_CTU.h"
#include "fwdmodel_dce_2CXM.h"

int main(int argc, char** argv) {

  //add the DCE models - these will autoregister at this point
  DCEFwdModel::NewInstance();
  DCE_ETM_FwdModel::NewInstance();
  DCE_CTU_FwdModel::NewInstance();
  DCE_2CXM_FwdModel::NewInstance();
  
  return execute(argc, argv);

}