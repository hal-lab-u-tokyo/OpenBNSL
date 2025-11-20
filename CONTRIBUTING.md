We are happy to accept contributions to the project in the form of pull requests.

# Issues
Please feel free to create a new issue for any bugs, questions etc.

---
# Scripts
```bash
pytest # run pytest to check all tests
./build_gtest.sh # build gtest (not used now)
./run_gtest.sh # run gtest (not used now)
./format.sh # format the code
./gen_docs.sh # build the documentation locally
```

--- 
# Pull Requests
1. Fork the repository
2. Create a new branch from the `develop` branch
3. Make your changes
4. Make sure all tests pass by running `pytest`
5. Format your code using the provided script: `./scripts/format.sh`
6. Check the documentation by running `./scripts/gen_docs.sh`
7. Push your changes to your fork
8. Create a pull request

---
# Branches
Branch naming examples:
- `master` - The main branch. This branch is always stable and contains the latest release.
- `develop` - The development branch. This branch contains the latest changes and is where new features are developed.
- `feature/*` - Feature branches. These branches are created from `develop` and are used to develop new features.

---
# Coding Style
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


---
# Interface of BNSL algorithm

The BNSL algorithm requires a `DataframeWrapper` class object as its first argument.  
The `DataframeWrapper` class provides the following interface:
``` C++
struct DataframeWrapper {
  DataframeWrapper(const py::object& dataframe); // Constructor from a Python dataframe
  ...
};
```

The output of the BNSL algorithm is a `PDAG` class object, which has the following interface:
``` C++
struct PDAG {
  std::size_t num_vars;
  std::vector<std::set<size_t>> parents; // parents[i]: parent set of node i
  ...
};
```

With this design, the learned PDAG can be easily converted to a pgmpy `PDAG` object as follows:

``` Python
import pandas as pd
from helpers.pgmpy_bridge import to_pgmpy

# Read a dataframe
dataframe = pd.read_csv("data.csv")
df_wrapper = openbnsl.base.DataframeWrapper(dataframe)

# Run the BNSL algorithm
learned_pdag_obnsl = openbnsl.structure_learning.hoge(df_wrapper, options)

# Convert to pgmpy PDAG
learned_pdag_pgmpy = to_pgmpy(learned_pdag_obnsl)
```
