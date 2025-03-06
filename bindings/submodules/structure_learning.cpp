#include <pybind11/pybind11.h>

#include "structure_learning/PC.h"
#include "structure_learning/RAI.h"
#include "structure_learning/cuPC.h"
#include "structure_learning/exhaustive_search.h"
#include "structure_learning/gpuPC.h"
#include "structure_learning/gpuPC2.h"
#include "structure_learning/gpuPC3.h"

void bind_structure_learning(py::module& m) {
  auto submodule =
      m.def_submodule("structure_learning", "Structure learning submodule");

  submodule.def("RAI", &RAI, "Run RAI algorithm");
  submodule.def("PC", &PC, "Run PC algorithm");
  submodule.def("exhaustive_search", &exhaustive_search,
                "Run exhaustive search algorithm", py::arg("df"),
                py::arg("score_type"), py::arg("max_parents"));
#ifdef USE_CUDA
  submodule.def("gpuPC", &cuda::gpuPC, "Run gpuPC algorithm");
  submodule.def("gpuPC2", &cuda2::gpuPC2, "Run gpuPC2 algorithm");
  submodule.def("gpuPC3", &cuda3::gpuPC3, "Run gpuPC3 algorithm");
  submodule.def("cuPC", &cuda_cupc::cuPC, "Run cuPC algorithm");
#endif
}
