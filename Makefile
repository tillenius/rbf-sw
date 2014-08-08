
FLAGS=-Wall -Ofast
INC=-I src/ -I src/mpi_superglue -I src/superglue
LIBS=-pthread

all: bin bin/swmpi bin/sw

bin:
	mkdir -p bin

bin/swmpi: bin
	mpic++ -DUSE_MPI -march=native ${FLAGS} src/standalone.cpp ${INC} ${LIBS} -o bin/swmpi

bin/sw: bin
	g++              -march=native ${FLAGS} src/standalone.cpp ${INC} ${LIBS} -o bin/sw

bin/swmpi.log: bin
	mpic++ -DUSE_MPI -DDEBUG_TRACE -march=native ${FLAGS} src/standalone.cpp ${INC} ${LIBS} -o bin/swmpi.log

bin/sw.log: bin
	g++              -DDEBUG_TRACE -march=native ${FLAGS} src/standalone.cpp ${INC} ${LIBS} -o bin/sw.log

.PHONY: all bin/swmpi bin/sw bin/swmpi.log bin/sw.log
