#include <pybind11/pybind11.h>

#include "score/score_type.h"

void bind_score(py::module& m) {
  auto submodule = m.def_submodule("score", "Score submodule");
  py::class_<BDeu>(submodule, "BDeu")
      .def(py::init<double>(), py::arg("ess") = 1.0)
      .def("__repr__", [](const BDeu& bdeu) {
        return "BDeu(" + std::to_string(bdeu.ess) + ")";
      });
  py::class_<K2>(submodule, "K2")
      .def(py::init<>())
      .def("__repr__", [](const K2& k2) { return "K2()"; });
}