CXX=g++-12
RM=rm -f
CFLAGS=-O3 -fopenmp

all: tsp-omp

tsp:
	$(CXX) tsp.cpp -o tsp $(CFLAGS)

tsp-omp:
	$(CXX) tsp-omp.cpp -o tsp-omp $(CFLAGS)

clean:
	$(RM) tsp-omp