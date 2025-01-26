#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "base/PDAG.h"
#include "base/dataframe_wrapper.h"
#include "score/score_type.h"
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
  // m.def("str2int_numpy", &str2int_numpy, "translate each string element of
  // numpy array into int element");
  m.def("PC", &PC, "Run PC algorithm");

  m.def("exhaustive_search", &exhaustive_search,
        "Run exhaustive search algorithm", py::arg("df"), py::arg("score_type"),
        py::arg("max_parents"));

  auto submodule_score = m.def_submodule("score", "Score submodule");
  py::class_<BDeu>(submodule_score, "BDeu")
      .def(py::init<double>(), py::arg("ess") = 1.0)
      .def("__repr__", [](const BDeu& bdeu) {
        return "BDeu(" + std::to_string(bdeu.ess) + ")";
      });
  py::class_<K2>(submodule_score, "K2")
      .def(py::init<>())
      .def("__repr__", [](const K2& k2) { return "K2()"; });

  auto submodule_base = m.def_submodule("base", "Base submodule");
  py::class_<DataframeWrapper>(submodule_base, "DataframeWrapper")
      .def(py::init<const py::object&>())
      .def("__repr__", [](const DataframeWrapper& df) {
        return "DataframeWrapper with " + std::to_string(df.num_of_datapoints) +
               " samples of " + std::to_string(df.num_of_vars) + " variables";
      });

#ifdef USE_CUDA
  m.attr("cuda_enabled") = true;
  m.def("matmul_cuda", &matmul_cuda, "Multiply two NumPy arrays with CUDA");
#endif
}
