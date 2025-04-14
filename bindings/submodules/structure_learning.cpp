#include <pybind11/pybind11.h>

// #include "structure_learning/PC_nishikori.h"
// #include "structure_learning/RAI_nishikori.h"
#include "structure_learning/PC.h"
#include "structure_learning/exhaustive_search.h"

void bind_structure_learning(py::module& m) {
  auto submodule =
      m.def_submodule("structure_learning", "Structure learning submodule");

  // submodule.def("PC_nishikori", &PC_nishikori, "Run PC Nishikori algorithm");
  // submodule.def("RAI_nishikori", &RAI_nishikori, "Run RAI Nishikori
  // algorithm");
  submodule.def("exhaustive_search", &exhaustive_search,
                "Run exhaustive search algorithm", py::arg("df"),
                py::arg("score_type"), py::arg("max_parents"));
  submodule.def("PC", &PC, "Run PC algorithm");
}
