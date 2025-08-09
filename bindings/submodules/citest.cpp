
#include <pybind11/pybind11.h>

#include "citest/citest_type.h"
#include "citest/citest.h"

void bind_citest(py::module& m) {
  auto submodule = m.def_submodule("citest", "CITests submodule");

  py::class_<ChiSquare>(submodule, "ChiSquare")
      .def(py::init<double>(), py::arg("level") = 0.05)
      .def("__repr__", [](const ChiSquare& chi_square) {
        return "ChiSquare(level=" + std::to_string(chi_square.level) + ")";
      });

  py::class_<GSquare>(submodule, "GSquare")
      .def(py::init<double>(), py::arg("level") = 0.05)
      .def("__repr__", [](const GSquare& g_square) {
        return "GSquare(level=" + std::to_string(g_square.level) + ")";
      });

  py::class_<OracleGraph>(submodule, "OracleGraph")
      .def(py::init<const PDAG&>(), py::arg("pdag"))
      .def("__repr__", [](const OracleGraph&) { return "<OracleGraph>"; });

  submodule.def("citest",
                &citest<true>,
                py::arg("x"),
                py::arg("y"),
                py::arg("sepset_candidate"),
                py::arg("ct"),
                py::arg("citest_type"));
}