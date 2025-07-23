#pragma once

#include "base/dataframe_wrapper.h"
#include "citest/citest_type.h"
#include "graph/PDAG.h"

/**
 * @brief A function to run the RAI algorithm
 *
 * @param df A DataframeWrapper object
 * @param max_parents The maximum number of parents for each node
 * @return A PDAG object
 */
PDAG RAI(const DataframeWrapper& df,
         const CITestType& ci_test_type,
         size_t max_cond_vars);