#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;


//void RAI(const py::array_t<std::string>& data, double ESS, py::array_t<bool>& g);
void RAI(const py::array_t<std::array<char, 8>>& data, double ESS, py::array_t<bool>& g);