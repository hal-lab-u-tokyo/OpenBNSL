
# Contributing 

We are happy to accept contributions to the project in the form of pull requests.

## Issues

Please feel free to create a new issue for any bugs, questions etc. 
It is very helpful if you gives us enough information to reproduce the problem. 
Github's guide on [about issues](https://guides.github.com/features/issues/) is also useful.

## Contributing Code

1. Fork the repository
2. Create a new branch from the `develop` branch
3. Make your changes
4. Push your changes to your fork
5. Create a pull request

# Scripts
```bash
pytest # run the frontend test
./build_test.sh # build the backend test
./run_test.sh # run the backend test
./format.sh # format the code
./gen_docs.sh # generate the documentation
```

## Branches

Branch naming examples:
- `master` - The main branch. This branch is always stable and contains the latest release.
- `develop` - The development branch. This branch contains the latest changes and is where new features are developed.
- `feature/*` - Feature branches. These branches are created from `develop` and are used to develop new features.
- `bugfix/*` - Bugfix branches. These branches are created from `develop` and are used to fix bugs.
- `hotfix/*` - Hotfix branches. These branches are created from `master` and are used to fix critical bugs.
- `docs/*` - Documentation branches. These branches are created from `develop` and are used to update documentation.
- `chore/*` - Chore branches. These branches are created from `develop` and are used for miscellaneous tasks.

## Directory Structure

```plaintext
.
├── .github/            # Github actions (build, test, format, etc.)
├── benchmarks/         # Benchmarking code
├── bindings/           # Python binding code (Pybind11) for the backend code
├── docs/               # Documentation (Doxygen, Sphinx)
├── examples/           # Jupyter notebook examples
├── external/           # External libraries (only pybind11 for now)
├── gtest/              # gtest test code for the backend code
├── images/             # Images for the README
├── include/            # Header files for the backend code
├── modules/            # Python modules
├── pytest/             # Pytest test code for the frontend code (test_*.py)
├── scripts/            # Scripts (format.sh, gen_docs.sh, etc.)
├── src/                # Source code for the backend code
├── .clang-format       # Clang format file
├── .gitignore          # Git ignore file
├── .gitmodules         # Git submodules file
├── CMakeLists.txt      # CMake build file (Compile backend code, bindings, and tests)
├── CONTRIBUTING.md     # Contributing guidelines
├── docker-compose.yml  # Docker compose file
├── Dockerfile          # Docker file
├── LICENSE             # License file
├── pyproject.toml      # Python project file
├── pytest.ini          # Pytest configuration file
├── README.md           # Project README
└── setup.py            # Python setup file
```

## Coding Style

For C++ code, we use clang-format with Google's C++ style guide [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
Before submitting a pull request, please run the following command to format your code:
```bash
./scripts/format.sh
```

| **Element** | **Naming Rule** | **Example** |
| --- | --- | --- |
| **File Name** | snake_case | `my_file.h`, `my_file.cpp` |
| **Directory Name** | snake_case | `src/`, `include/`, `test_utils/` |
| **Class** | PascalCase | `MyClass`, `MyStruct` |
| **Function** | snake_case | `my_function()` |
| **Variable** | snake_case | `my_variable` |
| **Constant** | UPPERCASE_SNAKE_CASE | `MY_CONSTANT` |
| **Namespace** | snake_case | `my_namespace` | `my_namespace::my_function()` |


## Data exchange between Python and C++
<!-- 
using pybind11 to create bindings for the C++ code
since expected data types is pandas DataFrame, 
we set dataframe wrapper class
all bnsl functions has df_wrapper as input
expected output is pgmpy graph object
the detail is TBD
 -->


## Data exchange between R and C++
<!--
using Rcpp to create bindings for the C++ code
since expected data types is data.frame,
we set dataframe wrapper class
all bnsl functions has df_wrapper as input
expected output is bnlearn graph object
the detail is TBD
-->

## Discussion
<!-- TBD --> 
