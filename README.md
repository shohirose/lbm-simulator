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
$ ./build/simulator poiseuille -f data/poiseuille.json
$ ./build/simulator cavity -f data/cavity.json
```