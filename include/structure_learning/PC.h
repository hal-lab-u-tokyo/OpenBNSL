#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "base/PDAG.h"
namespace py = pybind11;

PDAGwithAdjMat PC(const py::array_t<int> &data);
