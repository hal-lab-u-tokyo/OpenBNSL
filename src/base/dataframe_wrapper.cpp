#include "base/dataframe_wrapper.h"

#include <set>

DataframeWrapper::DataframeWrapper(const py::object& dataframe) {
  // get the column names and assign indices in lexicographical order
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
  py::array numpy_array = dataframe.attr("values");
  auto dataset_buf = numpy_array.request();
  if (dataset_buf.ndim != 2)
    throw std::invalid_argument("Input must be a 2D array");
  num_of_vars = dataset_buf.shape[1];
  num_of_datapoints = dataset_buf.shape[0];
  const py::object* dataset_ptr = static_cast<py::object*>(dataset_buf.ptr);

  // get the unique values for each column and assign indices in lexicographical
  // order
  val_str2idx.resize(num_of_vars);
  val_idx2str.resize(num_of_vars);
  for (size_t i = 0; i < num_of_vars; i++) {
    std::set<std::string> unique_values;
    for (size_t j = 0; j < num_of_datapoints; j++) {
      std::string value_str;
      try {
        value_str = dataset_ptr[i * num_of_datapoints + j].cast<std::string>();
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

  // manage the num_of_values
  num_of_values.resize(num_of_vars);
  for (size_t i = 0; i < num_of_vars; i++) {
    num_of_values[i] = val_idx2str[i].size();
  }

  // manage the data_column_major
  data_column_major.resize(num_of_vars);
  for (size_t i = 0; i < num_of_vars; i++) {
    data_column_major[i].resize(num_of_datapoints);
    for (size_t j = 0; j < num_of_datapoints; j++) {
      auto value = dataset_ptr[i * num_of_datapoints + j].cast<std::string>();
      data_column_major[i][j] = val_str2idx[i][value];
    }
  }

  // manage the data_row_major
  data_row_major.resize(num_of_datapoints);
  for (size_t i = 0; i < num_of_datapoints; i++) {
    data_row_major[i].resize(num_of_vars);
    for (size_t j = 0; j < num_of_vars; j++) {
      data_row_major[i][j] = data_column_major[j][i];
    }
  }

  // // debug print
  // std::cout << "num_of_vars: " << num_of_vars << std::endl;
  // std::cout << "num_of_datapoints: " << num_of_datapoints << std::endl;
  // for (size_t i = 0; i < num_of_vars; i++) {
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
