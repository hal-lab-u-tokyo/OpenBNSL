*Note: This repository is currently under active development. Please be aware that functionalities and features are still being finalized.*

1. [More About OpenBNSL](#more-about-openbnsl)
2. [Set Up](#set-up)
    1. [Using Docker](#using-docker)
        1. [Optional: Nvidia GPU Support](#optional-nvidia-gpu-support)
        2. [Optional: Gurobi License](#optional-gurobi-license)
    2. [Local Installation](#local-installation)
3. [Build and Install](#build-and-install)
4. [Features](#features)
5. [Tutorial and Examples](#tutorial-and-examples)
6. [Contributing](#contributing)
7. [License](#license)
8. [Acknowledgments](#acknowledgments)

# More About OpenBNSL
OpenBNSL is a framework for Benchmarking Bayesian Network Structure Learning (BNSL) methods.

⚠️ Under construction ⚠️

# Set Up

## (Recommended) Using Docker

```bash
git clone --recurse-submodules
cd openbnsl
docker compose build
docker compose up
```

### Optional: Nvidia GPU Support
If you have an Nvidia GPU, you can use the following command to build the image with GPU support.
```bash
docker compose build --build-arg BASE_IMAGE=nvidia/cuda:12.6.2-devel-ubuntu22.04
```

### Optional: Gurobi License 
If you want to use gurobi, you need to get a license file (gurobi.lic) and put it in the same directory as the Dockerfile.
Gurobi provides free licenses for academic use.
```bash
docker compose build --build-arg INSTALL_GUROBI=true
```

## (Not Recommended) Local Installation
⚠️ Under construction ⚠️

## Build and Install
```bash
pip install . # build and install the package
python3 setup.py build_ext --inplace # build in place
```

# Features
- Structure Learning
    - [x] Exhaustive Search
    - [x] Peter-Clark algorithm (PC)
    - [x] Recursive Autonomy Identification (RAI)
- Evaluation Metrics
    - [x] Structural Hamming Distance
    - [ ] Marginal Likelihood
        - [x] Bayesian Dirichlet equivalent uniform (BDeu) 
        - []
    - [ ] Inference Accuracy
- Benchmarking 
    - [ ] ⚠️ Under construction ⚠️

# Tutorial and Examples
⚠️ Under construction ⚠️

## Contributing
Please feel free to create a new issue for any bugs, questions etc. 
If you want to contribute code, please follow the instructions in [CONTRIBUTING](CONTRIBUTING.md)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments
This work was supported by
JSPS KAKENHI, Grant Number 24KJ0578,
JST CREST, Grant Number JPMJCR21D2, and
JST SPRING, Grant Number JPMJSP2108.