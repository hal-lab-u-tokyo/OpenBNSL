#pragma once

#include "base/dataframe_wrapper.h"
#include "citest/citest_type.h"
#include "graph/pdag.h"

/**
 * @ingroup structure_learning
 * @brief A function to run the RAI (Recursive Autonomy Identification)
 * algorithm.
 * @param df A DataframeWrapper object
 * @param ci_test_type The type of conditional independence test to use
 * @param max_cond_vars The maximum number of conditional variables to consider
 * @return A PDAG object
 */
PDAG rai(const DataframeWrapper& df,
         const CITestType& ci_test_type,
         size_t max_cond_vars,
         bool apply_meek_r4 = false);
