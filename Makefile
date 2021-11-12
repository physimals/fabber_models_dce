include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_dce
LIBS = -lfsl-fabberexec -lfsl-fabbercore -lfsl-newimage \
       -lfsl-miscmaths -lfsl-utils -lfsl-cprob \
       -lfsl-NewNifti -lfsl-znz -ldl
XFILES = fabber_dce
SOFILES = libfsl-fabbermodels_dce.so

# Forward models
OBJS =  fwdmodel_dce.o fwdmodel_dce_CTU.o fwdmodel_dce_2CXM.o fwdmodel_dce_AATH.o fwdmodel_dce_tofts.o
# Removed fwdmodel_dce_LLS.o fwdmodel_dce_Patlak.o fwdmodel_dce_ETM.o fwdmodel_dce_ETM_LLS.o fwdmodel_dce_CTU_LLS.o fwdmodel_dce_2CXM_LLS.o

# For debugging:
#OPTFLAGS = -ggdb

# Pass Git revision details
GIT_SHA1:=$(shell git describe --match=NeVeRmAtCh --always --abbrev=40 --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#

all: ${XFILES} ${SOFILES}

# models in a library
libfsl-fabbermodels_dce.so : ${OBJS}
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ ${LDFLAGS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_dce : fabber_client.o | libfsl-fabbermodels_dce.so
	${CXX} ${CXXFLAGS} -o $@ $< -lfsl-fabbermodels_dce ${LDFLAGS}
