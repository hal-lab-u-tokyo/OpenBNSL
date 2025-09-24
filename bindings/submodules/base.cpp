#include <pybind11/pybind11.h>

#include "base/contingency_table.h"
#include "base/dataframe_wrapper.h"
#include "graph/pdag.h"

void bind_base(py::module& m) {
  auto submodule = m.def_submodule("base", "Base submodule");

  py::class_<DataframeWrapper>(submodule, "DataframeWrapper")
      .def_readonly("num_vars", &DataframeWrapper::num_vars)
      .def_readonly("num_datapoints", &DataframeWrapper::num_datapoints)
      .def_readonly("col_idx2str", &DataframeWrapper::col_idx2str)
      .def_readonly("col_str2idx", &DataframeWrapper::col_str2idx)
      .def_readonly("val_idx2str", &DataframeWrapper::val_idx2str)
      .def_readonly("val_str2idx", &DataframeWrapper::val_str2idx)
      .def_readonly("num_of_values", &DataframeWrapper::num_of_values)
      .def_readonly("data_column_major", &DataframeWrapper::data_column_major)
      .def_readonly("data_row_major", &DataframeWrapper::data_row_major)
      .def(py::init<const py::object&>())
      .def("__repr__", [](const DataframeWrapper& df) {
        std::string repr = "DataframeWrapper(\n";
        repr +=
            "  num_datapoints: " + std::to_string(df.num_datapoints) + ",\n";
        repr += "  num_vars: " + std::to_string(df.num_vars) + ",\n";
        repr += "  columns: [";
        for (size_t i = 0; i < df.col_idx2str.size(); ++i) {
          repr += df.col_idx2str[i];
          if (i != df.col_idx2str.size() - 1) repr += ", ";
        }
        repr += "]\n)";
        return repr;
      });

  // only for testing purposes
  py::class_<ContingencyTable>(submodule, "ContingencyTable")
      .def(py::init<const std::vector<size_t>&, const DataframeWrapper&>(),
           py::arg("var_ids"),
           py::arg("df"))
      .def("marginalize_to",
           &ContingencyTable::marginalize_to,
           py::arg("var_ids_tgt"))
      .def_readonly("var_ids", &ContingencyTable::var_ids)
      .def_readonly("cardinalities", &ContingencyTable::cardinalities)
      .def_readonly("counts", &ContingencyTable::counts);

  py::class_<PDAG>(submodule, "PDAG")
      .def(py::init<size_t>())
      .def(py::init<const PDAG&>())
      .def_readonly("num_vars", &PDAG::num_vars)
      .def("has_edge", &PDAG::has_edge, py::arg("from"), py::arg("to"))
      .def("add_edge", &PDAG::add_edge, py::arg("from"), py::arg("to"))
      .def("remove_edge", &PDAG::remove_edge, py::arg("from"), py::arg("to"))
      .def("score", &PDAG::score, py::arg("df"), py::arg("score_type"))
      .def("__repr__", [](const PDAG& g) {
        std::string repr = "PDAG(\n";
        repr += "  num_vars: " + std::to_string(g.num_vars) + ",\n";
        repr += "  edges: [";
        bool first = true;
        for (size_t v = 0; v < g.num_vars; ++v) {
          for (auto u : g.parents[v]) {
            if (!first) repr += ", ";
            repr += std::to_string(u) + " -> " + std::to_string(v);
            first = false;
          }
        }
        repr += "]\n";
        repr += ")";
        return repr;
      });
}
