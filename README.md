# Overview

This is a collection of simulators using Lattice Boltzman Method. Currently, the following simulators are implemented:

- Poiseuille flow

# How to Build

CMake is required to build the programs. Please run the following command under the root directory:

```terminal
$ cmake -B build -S .
$ cmake --build build
```

# How to Run

To run `poiseuille-flow-2d`, type the following command under the project root directory:

```bash
$ build/poiseuille-flow-2d -f data/poiseuille-fow-2d.json
```