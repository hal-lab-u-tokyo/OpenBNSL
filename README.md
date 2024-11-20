*Note: This repository is a mock-up, and the content within is currently under active development. Please be aware that functionalities and features are still being finalized.*

# OpenBN: Open Bayesian Networks Library
- OpenBN is a Python library for structure learning of Bayesian Networks.
- The backend is written in C/C++, OpenMP, CUDA, (and HLS for FPGA in the future) to provide high performance and scalability.
- The frontend is written in Python, user-friendly, and easy to use, making use of ample libraries for data processing and visualization.
- the joint is done with [pybind11](https://github.com/pybind/pybind11), which allows for seamless integration between backend and frontend.
   - So, clone the repository with `-recurse-submodules` to include the submodules.
- The library is still under development and many features are still missing. 
- The goal is to provide a unified evaluation framework for structure learning of Bayesian Networks.

## Table of Contents
- [Installation](#installation)
- [Planed Features](#planed-features)
- [Set Up Docker Environment (Recommended)](#set-up-docker-environment)
- [Set Up Bare Metal Environment (Not Recommended)](#set-up-bare-metal-environment)
- [Build the project](#build-the-project)
- [Run Tests](#run-tests)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Planed Features
- [ ] Create Bayesian Networks manually
- [ ] Structure Learning from data
    - [ ] RAI algorithm
    - [ ] PC algorithm
    - [ ] Greedy Hill Climbing
    - [ ] Tabu Search
    - [ ] Simulated Annealing
    - [ ] Integer Programming
- [ ] Parameter Learning from data
- [ ] Inference with Bayesian Networks
- [ ] Visualization of Bayesian Networks
- [ ] Export and Import of Bayesian Networks
- [ ] Integration with other libraries


## Set Up Docker Environment (Recommended)

### Optional: Gurobi License
If you want to use gurobi, you need to get a license file (gurobi.lic) and put it in the same directory as the Dockerfile.
Gurobi provides free licenses for academic use.

### Build Docker Image
```bash
docker-compose build
```
or with gurobi
```bash
docker-compose build --build-arg INSTALL_GUROBI=true
```

### Run Docker Container
```bash
docker-compose up
```

## Set Up Bare Metal Environment (Not Recommended)

Please make sure you have the following dependencies installed on your system:
- C++17
- CMake 3.22
- Python 3.10
- Pandas 2.2.3
- NetworkX 3.4.2 
- Matplotlib 3.9.2
- Jupyter Notebook 7.2.2
- pyLPsolve 2.9.0
- (optional) Gurobi 11.0.0

## Build the project
```bash
./do_build.sh
```
compile shared library, python bindings, and backend tests.

## Run Tests
```bash
./do_test.sh
```
backend tests (Google Test) and frontend tests (pytest) will be run.

## Examples
<!-- TDB -->

## Contributing
<!-- TDB -->

## License
<!-- TDB -->

## Acknowledgments
This work was supported by
JSPS KAKENHI, Grant Number 24KJ0578,
JST CREST, Grant Number JPMJCR21D2, and
JST SPRING, Grant Number JPMJSP2108.