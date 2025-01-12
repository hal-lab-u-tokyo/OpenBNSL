#pragma once
#include <pybind11/pybind11.h>

#include "base/dataframe_wrapper.h"
namespace py = pybind11;

/**
 * @brief A function to run the exhaustive search with dynamic programming
 *
 * @param df A DataframeWrapper object
 * @param max_parents The maximum number of parents for each node
 * @return A NumPy array of boolean values
 *
 * @details
 * This function runs the exhaustive search with dynamic programming.
 * The time complexity of this algorithm is O(n^2 * 2^n) and the space
 * complexity is O(n * 2^n) where n is the number of nodes. When n is large,
 * this algorithm is not recommended to use.
 */
py::array_t<bool> exhaustive_search(const DataframeWrapper& df,
                                    int max_parents);

double score(const std::vector<int>& bucket_size);