#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

namespace cuda {
py::array_t<bool> gpuPC(py::array_t<uint8_t> data, py::array_t<int> n_states,
                        py::array_t<bool> true_model);
}  // namespace cuda