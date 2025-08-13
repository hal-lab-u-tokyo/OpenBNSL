# Scripts

- [`build_gtest.sh`]: Configures and builds the C++ project with `BUILD_TESTS=ON`, producing the Google Test executable (`openbnsl_test`) in the `./build/` directory.
- [`run_gtest.sh`]: Runs the Google Test executable (`openbnsl_test`) from the `./build/` directory if it exists.
- [`format.sh`]: Formats C/C++ and Python source files.
  - C/C++: Uses `clang-format` with rules from `.clang-format`.
  - Python: Uses `black` with configuration from `pyproject.toml`.
  Supports `--check` mode for verifying formatting without modifying files.
  Target directories:
    - C/C++: `src/`, `include/`, `gtest/`, `bindings/`  
    - Python: `benchmarks/`, `pytest/`, `helpers/`
- [`gen_docs.sh`]: Generates the project documentation: