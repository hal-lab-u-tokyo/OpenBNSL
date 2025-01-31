#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <map>
#include <string>
#include <vector>
namespace py = pybind11;

using dtype = uint8_t;

/**
 * @class DataframeWrapper
 * @brief A wrapper class for a pandas DataFrame.
 *
 * This class is used to store the data from a pandas DataFrame.
 * It stores the data in two formats: column-major and row-major.
 * It also stores the mapping between the column names and the column indices,
 * and the mapping between the column values and the column value indices.
 * using uint8_t to minimize data size and improve cache hit rate
 * (uint8_t is still redundant)
 */
class DataframeWrapper {
 public:
  DataframeWrapper(const py::object& dataframe);
  size_t num_of_vars;
  size_t num_of_datapoints;

  std::vector<std::string> col_idx2str;
  std::map<std::string, size_t> col_str2idx;
  std::vector<std::vector<std::string>> val_idx2str;
  std::vector<std::map<std::string, dtype>> val_str2idx;

  std::vector<size_t> num_of_values;
  std::vector<std::vector<dtype>>
      data_col_idx;  // [num_of_column][num_of_datapoint]
  std::vector<std::vector<dtype>>
      data_idx_col;  // [num_of_datapoint][num_of_column]
};
