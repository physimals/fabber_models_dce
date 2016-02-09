cmake_minimum_required(VERSION 3.2)
project(fabber_models_dce)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    fabber_client.cc
    fwdmodel_dce.cc
    fwdmodel_dce.h
    fwdmodel_dce_2CXM.cc
    fwdmodel_dce_2CXM.h
    fwdmodel_dce_2CXM_LLS.cc
    fwdmodel_dce_2CXM_LLS.h
    fwdmodel_dce_AATH.cc
    fwdmodel_dce_AATH.h
    fwdmodel_dce_CTU.cc
    fwdmodel_dce_CTU.h
    fwdmodel_dce_CTU_LLS.cc
    fwdmodel_dce_CTU_LLS.h
    fwdmodel_dce_ETM.cc
    fwdmodel_dce_ETM.h
    fwdmodel_dce_ETM_LLS.cc
    fwdmodel_dce_ETM_LLS.h
    fwdmodel_dce_LLS.cc
    fwdmodel_dce_LLS.h
    fwdmodel_dce_Patlak.cc
    fwdmodel_dce_Patlak.h)

add_executable(fabber_models_dce ${SOURCE_FILES})