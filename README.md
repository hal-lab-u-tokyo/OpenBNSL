*Note: This repository is a mock-up, and the content within is currently under active development. Please be aware that functionalities and features are still being finalized.*

# OpenBN: Open Bayesian Networks Library

## Planed Features
- [ ] Create Bayesian Networks manually
- [ ] Structure Learning from data
    - [ ] RAI algorithm
    - [ ] PC algorithm
    - [ ] Greedy Hill Climbing
    - [ ] Tabu Search
    - [ ] Simulated Annealing
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
- C++ compiler with C++?? support
- CMake version 3.?? or higher
- Python version 3.?? or higher
- Pandas version ?? or higher
- NetworkX version ?? or higher
- Matplotlib version ?? or higher
- Jupyter Notebook version ?? or higher
- pyLPsolve version ?? or higher
- (optional) Gurobi version ?? or higher


### Build and Install OpenBN
1. Navigate to the project directory.
2. Create a build directory:
   ```bash
   mkdir build
   cd build
   ```
3. Configure the project with CMake:
   ```bash
   cmake ..
   ```
4. Build the project:
   ```bash
   make
   ```

### Run Tests
<!-- TDB -->

### Examples
<!-- TDB -->
