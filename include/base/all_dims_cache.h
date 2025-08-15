#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <boost/multiprecision/cpp_int.hpp>
#include <memory>
#include <vector>

#include "base/dataframe_wrapper.h"
namespace py = pybind11;
namespace mp = boost::multiprecision;

struct Node {
 public:
  std::vector<int> freq_tbl;
  std::vector<std::unique_ptr<Node>> children;
};

/**
 * @ingroup base
 * @struct AllDimsCache
 * @brief Represents a cache for frequency tables of all dimensions given df.
 */
struct AllDimsCache {
  const DataframeWrapper& df;
  size_t max_depth;

  std::unique_ptr<Node> root_node;
  void branch(Node* parent,
              const size_t curr_depth,
              const int rightmost_var,
              const std::vector<int>& freq_tbl,
              const std::vector<int>& indices);

  AllDimsCache(const DataframeWrapper& df, int max_depth);
  const std::vector<int>& get_freq_tbl(const std::vector<int>& variables);
};