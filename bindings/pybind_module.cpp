#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "base/PDAG.h"
#include "structure_learning/PC.h"
#include "structure_learning/RAI.h"
#include "template.h"
namespace py = pybind11;

PYBIND11_MODULE(openbnsl, m) {
  m.doc() = "OpenBNSL Python bindings";
  m.attr("cuda_enabled") = false;

  m.def("matmul_naive", &matmul_naive, "Multiply two NumPy arrays");
  m.def("matmul_openmp", &matmul_openmp,
        "Multiply two NumPy arrays with OpenMP");
  m.def("RAI", &RAI, "Run RAI algorithm");
  //m.def("str2int_numpy", &str2int_numpy, "translate each string element of numpy array into int element");
  m.def("PC", &PC, "Run PC algorithm");
#ifdef USE_CUDA
  m.attr("cuda_enabled") = true;
  m.def("matmul_cuda", &matmul_cuda, "Multiply two NumPy arrays with CUDA");
#endif
}
