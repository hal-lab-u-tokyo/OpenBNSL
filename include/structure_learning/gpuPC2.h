#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

namespace cuda2 {
py::array_t<bool> gpuPC2(py::array_t<uint8_t> data, py::array_t<int> n_states);
}  // namespace cuda2