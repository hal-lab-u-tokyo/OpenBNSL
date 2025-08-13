# Scripts

- [`./build_gtest.sh`](./build_gtest.sh): Builds the the Google Test executable (`openbnsl_test`) in the `./build/` directory.
- [`./run_gtest.sh`](./run_gtest.sh): Runs the Google Test executable (`openbnsl_test`) from the `./build/` directory if it exists.
- [`./format.sh`](./format.sh): Formats C/C++ and Python source files using `clang-format` and `black` respectively (available `--check` mode).
    - C/C++: `src/`, `include/`, `gtest/`, `bindings/`
    - Python: `benchmarks/`, `pytest/`, `helpers/`
- [`./gen_docs.sh`](./gen_docs.sh): Generates the project documentation to `../docs/sphinx/build/html` using Doxygen and Sphinx.