include ${FSLCONFDIR}/default.mk

PROJNAME=eddy

USRCXXFLAGS=${PARALLELFLAGS}
USRINCFLAGS=-I${INC_NEWMAT} -I${INC_BOOST} 
USRLDFLAGS=-L${LIB_NEWMAT} -rdynamic ${PARALLELFLAGS}

OBJS=EddyCommandLineOptions.o ECModels.o ECScanClasses.o EddyUtils.o EddyHelperClasses.o DiffusionGP.o b0Predictor.o ../topup/topup_file_io.o ../topup/displacement_vector.o
LIBS=-lwarpfns -lmeshclass -lbasisfield -lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lutils -lnewmat -lm -lz

XFILES=eddy test_bench

all: eddy

eddy: eddy.o ${OBJS}
	${CXX} ${LDFLAGS} -o $@ eddy.o ${OBJS} ${LIBS}

test_bench: test_bench.o ${OBJS}
	${CXX} ${LDFLAGS} -o $@ test_bench.o ${OBJS} ${LIBS}

