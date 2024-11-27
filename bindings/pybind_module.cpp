#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "template.h"
#include "base/DAG.h"
namespace py = pybind11;

PYBIND11_MODULE(openbn, m) {
    m.doc() = "OpenBN Python bindings";
    
    m.def("matmul_naive", &matmul_naive, "Multiply two NumPy arrays");
    m.def("matmul_openmp", &matmul_openmp, "Multiply two NumPy arrays with OpenMP");

    py::class_<DAG>(m, "DAG")
        .def(py::init<size_t>())
        .def("get_adj_list", &DAG::get_adj_list);
    
    
}
