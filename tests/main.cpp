#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

#include <gtest/gtest.h>

namespace py = pybind11;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    py::scoped_interpreter guard{};
    return RUN_ALL_TESTS();
}