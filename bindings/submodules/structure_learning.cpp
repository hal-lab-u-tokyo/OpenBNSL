#include <pybind11/pybind11.h>

#include "structure_learning/exhaustive_search.h"
#include "structure_learning/pc.h"
#include "structure_learning/rai.h"
#include "structure_learning/simulated_annealing.h"

void bind_structure_learning(py::module& m) {
  auto submodule =
      m.def_submodule("structure_learning", "Structure learning submodule");

  submodule.def("exhaustive_search",
                &exhaustive_search,
                "Run exhaustive search algorithm",
                py::arg("df"),
                py::arg("score_type"),
                py::arg("max_parents"),
                py::arg("is_deterministic") = false);
  submodule.def("simulated_annealing",
                &simulated_annealing,
                "Run simulated annealing algorithm",
                py::arg("df"),
                py::arg("score_type"),
                py::arg("max_parents"),
                py::arg("max_iters") = 10000,
                py::arg("init_temp") = 1.0,
                py::arg("cooling_rate") = 0.9995,
                py::arg("is_deterministic") = false,
                py::arg("seed") = 0,
                py::arg("num_chains") = 0);
  submodule.def("pc",
                &pc,
                "Run PC algorithm",
                py::arg("df"),
                py::arg("ci_test_type"),
                py::arg("max_cond_vars"),
                py::arg("stable") = true);
  submodule.def("rai",
                &rai,
                "Run RAI algorithm",
                py::arg("df"),
                py::arg("ci_test_type"),
                py::arg("max_cond_vars"));
}
