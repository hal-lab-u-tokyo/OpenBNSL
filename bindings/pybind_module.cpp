#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "base/PDAG.h"
#include "template.h"
namespace py = pybind11;

PYBIND11_MODULE(openbnsl, m) {
  m.doc() = "OpenBNSL Python bindings";
  m.attr("cuda_enabled") = false;

  m.def("matmul_naive", &matmul_naive, "Multiply two NumPy arrays");
  m.def("matmul_openmp", &matmul_openmp,
        "Multiply two NumPy arrays with OpenMP");
#ifdef USE_CUDA
  m.attr("cuda_enabled") = true;
  m.def("matmul_cuda", &matmul_cuda, "Multiply two NumPy arrays with CUDA");
#endif
}
