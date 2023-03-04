CXX=g++-12
RM=rm -f
CFLAGS=-O3 -fopenmp -ggdb3

all: tsp

tsp:
	$(CXX) tsp.cpp -o tsp $(CFLAGS)

clean:
	$(RM) tsp

test: clean tsp
	python3 testing.py