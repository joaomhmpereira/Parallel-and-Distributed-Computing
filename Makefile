CXX=g++
RM=rm -f
CFLAGS=-O3 -fopenmp

all: tsp

tsp:
	$(CXX) tsp.cpp -o tsp $(CFLAGS)

clean:
	$(RM) tsp