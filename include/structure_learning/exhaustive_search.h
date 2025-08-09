#pragma once
#include <pybind11/pybind11.h>

#include "base/dataframe_wrapper.h"
#include "graph/pdag.h"
#include "score/score_type.h"
namespace py = pybind11;

/**
 * @brief A function to run the exhaustive search with dynamic programming
 *
 * @param df A DataframeWrapper object
 * @param max_parents The maximum number of parents for each node
 * @return A PDAG object representing the learned structure
 */

PDAG exhaustive_search(const DataframeWrapper& df,
                       const ScoreType& score_type,
                       size_t max_parents,
                       bool is_deterministic);