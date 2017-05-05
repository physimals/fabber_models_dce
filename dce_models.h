#pragma once

#ifdef _WIN32
#ifdef fabber_dce_EXPORTS
#define FABBER_DCE_API __declspec(dllexport)
#else
#define FABBER_DCE_API __declspec(dllimport)
#endif
#define CALL __stdcall
#else
#define FABBER_DCE_API
#define CALL
#endif

#include "fabber_core/fwdmodel.h"

extern "C"
{
FABBER_DCE_API int CALL get_num_models();
FABBER_DCE_API const char * CALL get_model_name(int index);
FABBER_DCE_API NewInstanceFptr CALL get_new_instance_func(const char *name);
}
