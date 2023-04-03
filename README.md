# Overview

This is a collection of simulators using Lattice Boltzman Method. Currently, the following simulators are implemented:

- 2-D Poiseuille flow
- 2-D cavity flow

# How to Build

CMake is required to build programs. Please run the following commands under the root directory:

```terminal
$ cmake -B build -S .
$ cmake --build build
```

# How to Run

Please type the following commands under the project root directory to run LBM simulators:

```terminal
$ ./build/lbm-simulator poiseuille -f data/poiseuille.json
$ ./build/lbm-simulator cavity -f data/cavity.json
```