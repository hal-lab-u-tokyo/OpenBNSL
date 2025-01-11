#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

py::array_t<bool> RAI(py::array_t<int> data, py::array_t<int> n_states, float ESS, int parallel, int threshold_DP);