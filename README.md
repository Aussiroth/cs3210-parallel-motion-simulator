# Parallel Particle Simulation
**Au Liang Jun (A0173593W), Yan Hong Yao Alvin (A0159960R)**  
CS3210 Parallel Computing AY19/20 Sem 1


## Directory
```
│   README.md
│   a1_report.pdf
│   particle-simulator.cpp
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
### Compilation
`g++ particle-simulator.cpp -fopenmp -std=c++11`

### Execution
`./a.out < input/random1000.txt`
