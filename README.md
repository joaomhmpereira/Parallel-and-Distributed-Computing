# `Parallel and Distributed Computing` #

This repository contains work developed in the [Parallel and Distributed Computing](https://fenix.tecnico.ulisboa.pt/disciplinas/CPD2/2022-2023/2-semestre/pagina-inicial) masters course at IST.

More specifically, it contains three (one sequential and two parallel) implementations of the [Travelling Salesman Problem](https://en.wikipedia.org/wiki/Travelling_salesman_problem) using the [Branch and Bound](https://en.wikipedia.org/wiki/Branch_and_bound) algorithm. 

Regarding the parallel implementations, one was done using [OpenMP](https://en.wikipedia.org/wiki/OpenMP) and the other using [MPI (Message Passing Interface)](https://en.wikipedia.org/wiki/Message_Passing_Interface).

Grade: 18/20.

Author | Github
-------|-------
Leonor Barreiros      | https://github.com/leonormbarreiros
Catarina Bento        | https://github.com/catarinab
João Pereira          | https://github.com/joaomhmpereira

### `Instructions to run each implementation` ###

All implementations have their own directory and each already include a Makefile to compile the program.

Directories `pub-instances` and `output` contain the available tests (input files) and expected outputs respectively.

After compiling, use the following commands to run the programs:
- Serial
    - `./tsp <path-to-input-file> <max-cost>`
- OpenMP
    - `./tsp-omp <path-to-input-file> <max-cost>`
- MPI
    - `mpirun -n <number-of-processes> tsp-mpi <path-to-input-file> <max-cost>`
