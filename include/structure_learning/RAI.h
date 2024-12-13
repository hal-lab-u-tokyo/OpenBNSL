#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;

std::vector<std::vector<int>> RAI(std::vector<std::vector<std::string>>& data, float ESS);