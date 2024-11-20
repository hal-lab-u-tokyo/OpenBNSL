#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "template.h"
namespace py = pybind11;

PYBIND11_MODULE(openbn, m) {
    m.doc() = "OpenBN Python bindings";
    
    m.def("matmul_naive", &matmul_naive, "Multiply two NumPy arrays");
    m.def("matmul_openmp", &matmul_openmp, "Multiply two NumPy arrays with OpenMP");

}
