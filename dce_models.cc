/* dce_models.cc 

Michael Chappell - IBME & FMRIB Analysis Group

Copyright (C) 2010-2011 University of Oxford */

/* CCOPYRIGHT  */

#include "fabber_core/fwdmodel.h"

#include <algorithm>

#include "fwdmodel_dce_2CXM.h"
#include "fwdmodel_dce_AATH.h"
#include "fwdmodel_dce_CTU.h"
#include "fwdmodel_dce_tofts.h"

#include "dce_models.h"

extern "C" {

int CALL get_num_models()
{
    return 4;
}

const char *CALL get_model_name(int index)
{
    switch (index)
    {
    case 0:
        return "dce_tofts";
        break;
    case 1:
        return "dce_2CXM";
        break;
    case 2:
        return "dce_AATH";
        break;
    case 3:
        return "dce_CTU";
        break;
    default:
        return NULL;
    }
}

NewInstanceFptr CALL get_new_instance_func(const char *name)
{
    if (string(name) == "dce_tofts")
    {
        return DCEStdToftsFwdModel::NewInstance;
    }
    else if (string(name) == "dce_2CXM")
    {
        return DCE_2CXM_FwdModel::NewInstance;
    }
    else if (string(name) == "dce_AATH")
    {
        return DCE_AATH_FwdModel::NewInstance;
    }
    else if (string(name) == "dce_CTU")
    {
        return DCE_CTU_FwdModel::NewInstance;
    }
    else
    {
        return NULL;
    }
}
}
