include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_dce

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -I..
# uncomment below to activate motion correction
#USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -D__OXASL -D__FABBER_MOTION
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L../fabber_core

LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat -lfslio -lniftiio -lznz -lz -ldl
# uncomment below to activate motion correction
#LIBS = -lwarpfns -lutils -lnewimage -lmiscmaths -lprob -lnewmat -lfslio -lniftiio -lznz -lz

XFILES = fabber_dce

# Forward models
OBJS =  fwdmodel_dce.o fwdmodel_dce_LLS.o fwdmodel_dce_Patlak.o fwdmodel_dce_ETM.o fwdmodel_dce_ETM_LLS.o fwdmodel_dce_CTU.o fwdmodel_dce_CTU_LLS.o fwdmodel_dce_2CXM.o fwdmodel_dce_2CXM_LLS.o fwdmodel_dce_AATH.o fwdmodel_dce_2CXM_NLLS.o

# For debugging:
OPTFLAGS = -ggdb
#OPTFLAGS =

# Pass Git revision details
GIT_SHA1:=$(shell git describe --match=NeVeRmAtCh --always --abbrev=40 --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#

all:	${XFILES} libfabbermodels_dce.a

# models in a library
libfabbermodels_dce.a : ${OBJS}
	${AR} -r $@ ${OBJS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_dce : fabber_client.o ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${OBJS} -lfabbercore -lfabberexec ${LIBS}

# DO NOT DELETE
