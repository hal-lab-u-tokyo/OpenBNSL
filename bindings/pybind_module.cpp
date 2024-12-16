#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "base/PDAG.h"
#include "template.h"
#include "structure_learning/addadd.h"
#include "structure_learning/RAI.h"
namespace py = pybind11;

PYBIND11_MODULE(openbnsl, m) {
  m.doc() = "OpenBNSL Python bindings";

  m.def("matmul_naive", &matmul_naive, "Multiply two NumPy arrays");
  m.def("matmul_openmp", &matmul_openmp,
        "Multiply two NumPy arrays with OpenMP");
  m.def("RAI", &RAI, "Run RAI algorithm");
}
