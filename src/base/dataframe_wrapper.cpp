#include "base/dataframe_wrapper.h"

#include <set>

DataframeWrapper::DataframeWrapper(const py::object& dataframe) {
  // a helper function to normalize a cell value to string
  auto normalize = [&](const py::object& obj) -> std::string {
    if (obj.is_none()) return "<NA>";
    // future: handle other types
    return obj.cast<std::string>();
  };

  // get the column names
  if (!py::hasattr(dataframe, "columns"))
    throw std::invalid_argument("Input must be a pandas dataframe");
  auto columns = dataframe.attr("columns");
  for (auto column : columns) {
    std::string column_str;
    try {
      column_str = column.cast<std::string>();
    } catch (const std::exception& e) {
      throw std::invalid_argument("Failed to cast column name to string");
    }
    col_str2idx[column_str] = col_idx2str.size();
    col_idx2str.push_back(column_str);
  }

  // get the numpy array
  if (!py::hasattr(dataframe, "values"))
    throw std::invalid_argument("Input must be a pandas dataframe");
  py::array arr = dataframe.attr("values");
  auto buf = arr.request();
  if (buf.ndim != 2) throw std::invalid_argument("Input must be a 2D array");
  num_vars = static_cast<size_t>(buf.shape[1]);
  num_datapoints = static_cast<size_t>(buf.shape[0]);
  const py::object* ptr = static_cast<py::object*>(buf.ptr);

  // get the unique values for each column
  val_str2idx.resize(num_vars);
  val_idx2str.resize(num_vars);
#pragma omp parallel for
  for (size_t i = 0; i < num_vars; i++) {
    std::set<std::string> unique_values;
    for (size_t j = 0; j < num_datapoints; j++) {
      std::string value_str;
      try {
        value_str = normalize(ptr[i * num_datapoints + j]);
      } catch (const std::exception& e) {
        throw std::invalid_argument("Failed to cast value to string");
      }
      unique_values.insert(value_str);
    }
    for (const auto& value : unique_values) {
      val_str2idx[i][value] = val_idx2str[i].size();
      val_idx2str[i].push_back(value);
    }
  }

  // manage the num_values
  num_values.resize(num_vars);
  for (size_t i = 0; i < num_vars; i++) {
    num_values[i] = val_idx2str[i].size();
  }

  // manage the data_column_major
  data_column_major.resize(num_vars);
#pragma omp parallel for
  for (size_t i = 0; i < num_vars; i++) {
    data_column_major[i].resize(num_datapoints);
    for (size_t j = 0; j < num_datapoints; j++) {
      auto value_str = normalize(ptr[i * num_datapoints + j]);
      data_column_major[i][j] = val_str2idx[i][value_str];
    }
  }

  // manage the data_row_major
  data_row_major.resize(num_datapoints);
  for (size_t i = 0; i < num_datapoints; i++) {
    data_row_major[i].resize(num_vars);
    for (size_t j = 0; j < num_vars; j++) {
      data_row_major[i][j] = data_column_major[j][i];
    }
  }

  // // debug print
  // std::cout << "num_vars: " << num_vars << std::endl;
  // std::cout << "num_datapoints: " << num_datapoints << std::endl;
  // for (size_t i = 0; i < num_vars; i++) {
  //   std::cout
  //   << "column_strs[" << i << "]: "
  //   << col_idx2str[i]
  //   << ", values["
  //   << i << "]: {";
  //   for (size_t j = 0; j < val_idx2str[i].size(); j++) {
  //     std::cout << val_idx2str[i][j] << ", ";
  //   }
  //   std::cout
  //   << "}("
  //   << val_idx2str[i].size()
  //   << ")" << std::endl;
  // }
}
