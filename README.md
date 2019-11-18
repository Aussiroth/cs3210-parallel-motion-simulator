# Parallel Particle Simulation
**Au Liang Jun (A0173593W), Yan Hong Yao Alvin (A0159960R)**  
CS3210 Parallel Computing AY19/20 Sem 1

---
**NOTE**

Our OpenMP and CUDA implementations have been updated to fix a bug in Assignment 1 Part 1 and to follow the output format specified more closely.

---

## Directory
```
│   README.md
│   a2_report.pdf
|   cuda-simulator.cu
|   mpi-simulator.cpp
│   omp-simulator.cpp
│
└───input
│   │   random1000.txt
│   │   random2000.txt
|   |   random3000.txt
|   |   ...
│   │   random10000.txt
│   
└───output
    │   out1000.txt
    │   out2000.txt
    |   out3000.txt
    |   ...
    │   out10000.txt
```
## Instructions
### OpenMP
#### Compilation
`g++ omp-simulator.cpp -fopenmp -std=c++11 -o omp-simulator`

#### Execution
`./omp-simulator < input/random1000.txt`

### CUDA
#### Compilation
`nvcc cuda-simulator.cu -std=c++11 -o cuda-simulator`

#### Execution
`./cuda-simulator < input/random1000.txt`

### OpenMPI
#### Compilation
`mpic++ mpi-simulator.cpp -std=c++11 -o mpi-simulator`

#### Execution
`mpirun -np 2 ./mpi-simulator < input/random1000.txt`