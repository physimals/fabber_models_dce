/* dce_models.cc 

Michael Chappell - IBME & FMRIB Analysis Group

Copyright (C) 2010-2011 University of Oxford */

/* CCOPYRIGHT  */

#include "fabber_core/fwdmodel.h"

#include <algorithm>

#include "fwdmodel_dce.h"
#include "fwdmodel_dce_2CXM.h"
#include "fwdmodel_dce_2CXM_LLS.h"
#include "fwdmodel_dce_AATH.h"
#include "fwdmodel_dce_CTU.h"
#include "fwdmodel_dce_CTU_LLS.h"
#include "fwdmodel_dce_ETM.h"
#include "fwdmodel_dce_ETM_LLS.h"
#include "fwdmodel_dce_LLS.h"
#include "fwdmodel_dce_Patlak.h"
#include "fwdmodel_dce_tofts.h"

#include "dce_models.h"

extern "C" {

int CALL get_num_models()
{
    return 11;
}

const char * CALL get_model_name(int index)
{
    switch (index)
    {
    case 0:
        return "dce";
        break;
    case 1:
        return "dce_2CXM";
        break;
    case 2:
        return "dce_2CXM_LLS";
        break;
    case 3:
        return "dce_AATH";
        break;
    case 4:
        return "dce_CTU";
        break;
    case 5:
        return "dce_CTU_LLS";
        break;
    case 6:
        return "dce_ETM";
        break;
    case 7:
        return "dce_ETM_LLS";
        break;
    case 8:
        return "dce_LLS";
        break;
    case 9:
        return "dce_Patlak";
        break;
    case 10:
        return "dce_tofts";
        break;
    default:
        return NULL;
    }
}

NewInstanceFptr CALL get_new_instance_func(const char *name)
{
    if (string(name) == "dce")
    {
        return DCEToftsFwdModel::NewInstance;
    }
    else if (string(name) == "dce_2CXM")
    {
        return DCE_2CXM_FwdModel::NewInstance;
    }
    else if (string(name) == "dce_2CXM_LLS")
    {
        return DCE_2CXM_LLS_FwdModel::NewInstance;
    }
    else if (string(name) == "dce_AATH")
    {
        return DCE_AATH_FwdModel::NewInstance;
    }
    else if (string(name) == "dce_CTU")
    {
        return DCE_CTU_FwdModel::NewInstance;
    }
    else if (string(name) == "dce_CTU_LLS")
    {
        return DCE_CTU_LLS_FwdModel::NewInstance;
    }
    else if (string(name) == "dce_ETM")
    {
        return DCE_ETM_FwdModel::NewInstance;
    }
    else if (string(name) == "dce_ETM_LLS")
    {
        return DCE_ETM_LLS_FwdModel::NewInstance;
    }
    else if (string(name) == "dce_LLS")
    {
        return DCE_LLS_FwdModel::NewInstance;
    }
    else if (string(name) == "dce_Patlak")
    {
        return DCE_Patlak_FwdModel::NewInstance;
    } 
    else if (string(name) == "dce_tofts")
    {
        return DCEStdToftsFwdModel::NewInstance;
    }
    else
    {
        return NULL;
    }
}
}
