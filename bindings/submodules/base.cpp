#include <pybind11/pybind11.h>

#include "base/PDAG.h"
#include "base/PDAG2.h"
#include "base/contingency_table.h"
#include "base/dataframe_wrapper.h"

void bind_base(py::module& m) {
  auto submodule = m.def_submodule("base", "Base submodule");

  py::class_<DataframeWrapper>(submodule, "DataframeWrapper")
      .def_readonly("num_of_vars", &DataframeWrapper::num_of_vars)
      .def_readonly("num_of_datapoints", &DataframeWrapper::num_of_datapoints)
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
        repr += "  num_of_datapoints: " + std::to_string(df.num_of_datapoints) +
                ",\n";
        repr += "  num_of_vars: " + std::to_string(df.num_of_vars) + ",\n";
        repr += "  columns: [";
        for (size_t i = 0; i < df.col_idx2str.size(); ++i) {
          repr += df.col_idx2str[i];
          if (i != df.col_idx2str.size() - 1) repr += ", ";
        }
        repr += "]\n)";
        return repr;
      });

  py::class_<ContingencyTable>(submodule, "ContingencyTable")
      .def_readonly("vars", &ContingencyTable::vars)
      .def_readonly("cardinalities", &ContingencyTable::cardinalities)
      .def_readonly("counts", &ContingencyTable::counts);

  submodule.def("buildContingencyTable", &buildContingencyTable,
                py::arg("vars"), py::arg("df"));

  py::class_<PDAG>(submodule, "PDAG")
      .def(py::init<size_t>())
      .def(py::init<const PDAG&>())
      .def("has_edge", &PDAG::has_edge, py::arg("from"), py::arg("to"))
      .def("add_edge", &PDAG::add_edge, py::arg("from"), py::arg("to"))
      .def("remove_edge", &PDAG::remove_edge, py::arg("from"), py::arg("to"));
}
