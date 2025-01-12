#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "base/PDAG.h"
#include "base/dataframe_wrapper.h"
#include "structure_learning/PC.h"
#include "structure_learning/RAI.h"
#include "structure_learning/exhaustive_search.h"
#include "template.h"
namespace py = pybind11;

PYBIND11_MODULE(openbnsl, m) {
  m.doc() = "OpenBNSL Python bindings";
  m.attr("cuda_enabled") = false;

  m.def("matmul_naive", &matmul_naive, "Multiply two NumPy arrays");
  m.def("matmul_openmp", &matmul_openmp,
        "Multiply two NumPy arrays with OpenMP");
  m.def("RAI", &RAI, "Run RAI algorithm");
  m.def("PC", &PC, "Run PC algorithm");

  py::class_<DataframeWrapper>(m, "DataframeWrapper")
      .def(py::init<const py::object&>());
  m.def("exhaustive_search", &exhaustive_search,
        "Run exhaustive search algorithm", py::arg("df"),
        py::arg("max_parents"));

#ifdef USE_CUDA
  m.attr("cuda_enabled") = true;
  m.def("matmul_cuda", &matmul_cuda, "Multiply two NumPy arrays with CUDA");
#endif
}
