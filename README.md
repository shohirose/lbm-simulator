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

JSON files under the `data` directory are input files for each simulator.
An example of the input JSON file for the `poiseuille` simulator is:

```yaml
{
    // Number of grids
    "gridShape": [
        21, // X
        21  // Y
    ],
    // External force vector
    // Constant and uniform in the region
    "externalForce": [
        0.00001, // X
        0.0      // Y
    ],
    "relaxationTime": 0.56,
    // Criteria to end while-loops
    "errorLimit": 1e-10,
    "printFrequency": 5000,
    "maxIteration": 1000000,
    // Directory to output result data (ux.txt)
    "outputDirectory": "./result/poiseuille"
}
```

That for the `cavity` simulator is:

```yaml
{
    "gridShape": [
        52,
        52
    ],
    // Wall velocity of the top boundary
    "wallVelocity": 0.1,
    // Reynolds number, which is used to compute relaxation time
    "reynoldsNumber": 100,
    "errorLimit": 1e-10,
    "printFrequency": 10000,
    "maxIteration": 50000,
    // Directory to output result data (ux.txt, uy.txt)
    "outputDirectory": "./result/cavity"
}
```