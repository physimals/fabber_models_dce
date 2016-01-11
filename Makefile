include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_dce

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB}

LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat -lfslio -lniftiio -lznz -lz

XFILES = fabber

# Forward models
OBJS =  fwdmodel_dce.o fwdmodel_dce_ETM.o fwdmodel_dce_CTU.o fwdmodel_dce_2CXM.o

# For debugging:
OPTFLAGS = -ggdb
#OPTFLAGS =

#
# Build
#

all:	${XFILES} libfabbermodels_dce.a

# models in a library
libfabbermodels_dce.a : ${OBJS}
	${AR} -r $@ ${OBJS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber : fabber_client.o ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${OBJS} -lfabbercore ${LIBS}

# DO NOT DELETE
