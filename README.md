**Note: This repository is currently under active development.**

OpenBNSL is an open framework designed to enable fair and highly reproducible comparisons of Bayesian Network Structure Learning (BNSL) algorithms. This framework provides an environment for systematically and fairly comparing various BNSL techniques, supporting the further advancement of Bayesian Network research.

| Component              | Description                                               |
|------------------------|-----------------------------------------------------------|
| Core Library           | Fast BNSL in C/C++ with OpenMP & CUDA                     |
| Evaluation Suite       | Python scripts for evaluating BNSL algorithms             |
| Experiment Environment | Docker-based environment for reproducible experiments     |

![OpenBNSL Architecture](images/architecture.png)


1. [Set Up](#set-up)
    1. [Using Docker](#using-docker)
    2. [Using Bare Metal](#using-bare-metal)
2. [Build and Install](#build-and-install)
3. [Features](#features)
4. [Contributing](#contributing)
5. [License](#license)
6. [Acknowledgments](#acknowledgments)


---
# Set Up

## Using Docker
```bash
git clone --recurse-submodules git@github.com:hal-lab-u-tokyo/OpenBNSL.git # for pybind11 submodule
cd OpenBNSL
docker compose build \
  [--build-arg BASE_IMAGE=nvidia/cuda:12.6.2-devel-ubuntu22.04] \  # Nvidia GPU
  [--build-arg INSTALL_R=true] \                                   # R & bnlearn
  [--build-arg INSTALL_GUROBI=true]                                # Gurobi (requires license)
docker compose up
```

## Using Bare Metal
⚠️ Under construction ⚠️

---
# Build and Install
```bash
pip install . # build and install the package
python3 setup.py build_ext --inplace # build in place
```

---
# Features

- Core Library (C/C++ with OpenMP & CUDA)
    - Score-based Structure Learning
        - [x] Exhaustive Search
    - Constraint-based Structure Learning
        - [x] Peter-Clark algorithm (PC) 
        - [ ] Recursive Autonomy Identification (RAI)
    - local-to-global learning
        - [ ] 
- Evaluation Suite (Python)
    - Evaluation Metrics
        - [x] Structural Hamming Distance
        - [x] Marginal Likelihood
        - [ ] Inference accuracy
        - [ ] memory usage
        - [ ] runtime
    - Benchmarking Senarios
        - [ ] scalability for number of variables
        - [ ] scalability for number of samples
    - [ ] Generate Synthetic Data
- Experiment Environment (Docker)
    - [x] OpenMP support (default)
    - [x] Nvidia GPU and CUDA toolkit support
    - [x] Gurobi support
    - [ ] R and `bnlearn` support

--- 
# Contributing
Please feel free to create a new issue for any bugs, questions etc. 
If you want to contribute code, please follow the instructions in [CONTRIBUTING](CONTRIBUTING.md)

--- 
# License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
# Acknowledgments
This work was supported by
JSPS KAKENHI, Grant Number 24KJ0578,
JST CREST, Grant Number JPMJCR21D2, and
JST SPRING, Grant Number JPMJSP2108.