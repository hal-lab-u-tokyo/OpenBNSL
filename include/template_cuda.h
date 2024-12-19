#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;

py::array_t<double> matmul_cuda(py::array_t<double> A, py::array_t<double> B);
