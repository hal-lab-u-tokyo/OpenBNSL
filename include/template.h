#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

py::array_t<double> matmul_naive(py::array_t<double> A, py::array_t<double> B);
py::array_t<double> matmul_openmp(py::array_t<double> A, py::array_t<double> B);
