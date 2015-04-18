
FLAGS=-Wall -Ofast
INC=-I src/ -I mpi-superglue -I superglue
LIBS=-pthread

all: bin bin/swmpi2.dev bin/swmpi2.prod bin/swmpi2.log bin/sw2.dev bin/sw2.prod bin/sw2.log

bin:
	mkdir -p bin

bin/swmpi2.dev: bin
	mpic++ -DRBFSW_DEBUG -DUSE_MPI -march=native ${FLAGS} src/standalone.cpp ${INC} ${LIBS} -o $@

bin/swmpi2.prod: bin
	mpic++ -DUSE_MPI -march=native ${FLAGS} src/standalone.cpp ${INC} ${LIBS} -o $@

bin/swmpi2.log: bin
	mpic++ -DUSE_MPI -DDEBUG_TRACE -march=native ${FLAGS} src/standalone.cpp ${INC} ${LIBS} -o $@

bin/sw2.dev: bin
	g++    -DRBFSW_DEBUG           -march=native ${FLAGS} src/standalone.cpp ${INC} ${LIBS} -o $@

bin/sw2.prod: bin
	g++              -march=native ${FLAGS} src/standalone.cpp ${INC} ${LIBS} -o $@

bin/sw2.log: bin
	g++              -DDEBUG_TRACE -march=native ${FLAGS} src/standalone.cpp ${INC} ${LIBS} -o $@

.PHONY: all bin/swmpi2.dev bin/swmpi2.prod bin/swmpi2.log bin/sw2.dev bin/sw2.prod bin/sw2.log
