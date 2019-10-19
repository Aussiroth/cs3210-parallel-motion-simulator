# Parallel Particle Simulation
**Au Liang Jun (A0173593W), Yan Hong Yao Alvin (A0159960R)**  
CS3210 Parallel Computing AY19/20 Sem 1


## Directory
```
│   README.md
│   a1_report.pdf
│   omp-simulator.cpp
|   cuda-simulator.cu
│
└───input
│   │   random1000.txt
│   │   random2000.txt
|   |   random3000.txt
|   |   random4000.txt
│   
└───output
    │   out1000.txt
    │   out2000.txt
    |   out3000.txt
    |   out4000.txt
```
## Instructions
### OpenMP
#### Compilation
`g++ omp-simulator.cpp -fopenmp -std=c++11 -o omp-simulator`

#### Execution
`./omp-simulator < input/random1000.txt`

### CUDA
#### Compilation
`nvcc cuda-simulator.cu -o cuda-simulator`

#### Execution
`./cuda-simulator < input/random1000.txt`
