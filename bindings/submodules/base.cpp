#include <pybind11/pybind11.h>

#include "base/PDAG.h"
#include "base/PDAG2.h"
#include "base/dataframe_wrapper.h"
#include "base/CItest.h"

void bind_base(py::module& m) {
  auto submodule = m.def_submodule("base", "Base submodule");

  py::class_<DataframeWrapper>(submodule, "DataframeWrapper")
      .def(py::init<const py::object&>())
      .def("__repr__", [](const DataframeWrapper& df) {
        return "DataframeWrapper with " + std::to_string(df.num_of_datapoints) +
               " samples of " + std::to_string(df.num_of_vars) + " variables";
      });
}
