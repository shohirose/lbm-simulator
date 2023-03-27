# Overview

This is a collection of simulators using Lattice Boltzman Method. Currently, the following simulators are implemented:

- 2-D Poiseuille flow
- 2-D cavity flow

# How to Build

CMake is required to build the programs. Please run the following command under the root directory:

```terminal
$ cmake -B build -S .
$ cmake --build build
```

# How to Run

Type the following commands under the project root directory to run simulators:

```bash
$ ./build/poiseuille-flow-2d -f ./data/poiseuille-flow-2d.json -o ./result/poiseuille-flow-2d
$ ./build/cavity-flow-2d -f ./data/cavity-flow-2d.json -o ./result/cavity-flow-2d
```