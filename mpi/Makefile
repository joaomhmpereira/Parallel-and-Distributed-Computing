CXX=mpic++
RM=rm -f
CFLAGS=-O3 -fopenmp -ggdb3

all: tsp

tsp:
	$(CXX) tsp-mpi.cpp -o tsp-mpi $(CFLAGS)

clean:
	$(RM) tsp