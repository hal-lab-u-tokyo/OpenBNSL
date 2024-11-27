*Note: This repository is a mock-up, and the content within is currently under active development. Please be aware that functionalities and features are still being finalized.*

# OpenBNSL: Open BNSL Framework
Bayesian Network Structure Learning (BNSL) is the task of learning the graph structure of a Bayesian Network from data.
It is widely recognized as a challenging problem due to its computational complexity, and the difficulties associated with handling high-dimensional data.
To address these challenges, various methods have been proposed.

OpenBNSL is a unified, open-source, and comprehensive framework for evaluating BNSL methods.
Researchers can leverage this framework to demonstrate the advantages of their algorithms, while general users can utilize BNSL techniques without requiring specialized expertise.

To fulfill its purpose, the framework is designed with the following requirements in mind:
1. Performance: The provided algorithms should be fast enough for effective benchmarking.
2. Scalability: The framework should handle large datasets and a significant number of random variables (e.g., over 1,000 variables), except where limited by algorithmic constraints.
3. Reproducibility: The framework should facilitate the easy reproduction of experiments.

To meet these requirements, OpenBNSL is designed as follows:
1. Backend: Written in C/C++ with OpenMP for high performance and scalability, with plans for further acceleration in the future.
2. Frontend: Implemented in Python to provide a user-friendly interface and leverage powerful libraries for:
    - Data manipulation (`Pandas`)
    - Visualization (`Matplotlib`, `NetworkX`)
    - Interactive prototyping (`Jupyter Notebook`)
    - Utilizing learned models (`pgmpy`)
3. Integration: The backend and frontend are connected using pybind11, enabling seamless and high-performance communication between Python and C++.
4. Containerization: The project is containerized using Docker, ensuring a reproducible environment, objective benchmarking, and ease of contribution from the community.

This framework is offered as an open-source project, and we actively welcome contributions from the community. We hope OpenBNSL becomes a valuable tool for researchers and practitioners working in the field of Bayesian Networks.

Note: When cloning this repository, please make sure to clone the repository with `--recurse-submodules` to include the pybind11 submodule.

## Features
- Structure Learning from data
    - Score-based Learning
        - [ ] Hill Climbing
        - [ ] Tabu Search
        - [ ] Simulated Annealing
        - [ ] Integer Programming
    - Constraint-based Learning
        - [ ] PC algorithm
        - [ ] RAI algorithm
- Evaluation Metrics
    - Structural Likelihood
        - [ ] Akaike Information Criterion (AIC)
        - [ ] Bayesian Information Criterion (BIC)
        - [ ] Minimum Description Length (MDL)
        - [ ] Bayesian Dirichlet equivalent uniform (BDeu) 
    - [ ] Structural Hamming Distance
    - [ ] Inference Accuracy
- Benchmarking 

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

## Tutorial and Examples
TBD

## Contributing
Please feel free to create a new issue for any bugs, questions etc. 
It is very helpful if you gives us enough information to reproduce the problem. 
Github's guide on [about issues](https://guides.github.com/features/issues/) is also useful.

If you want to contribute code, please follow the instructions in [Contributing](CONTRIBUTING.md)

## Documentation
TDB

## License
TBD

## Acknowledgments
This work was supported by
JSPS KAKENHI, Grant Number 24KJ0578,
JST CREST, Grant Number JPMJCR21D2, and
JST SPRING, Grant Number JPMJSP2108.