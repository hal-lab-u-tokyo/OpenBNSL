#pragma once

#include "base/PDAG.h"
#include "base/dataframe_wrapper.h"
#include "citest/citest_type.h"

/**
 * @brief A function to run the PC2 algorithm
 *
 * @param df A DataframeWrapper object
 * @param max_parents The maximum number of parents for each node
 * @return A PDAG object
 */
PDAG PC2(const DataframeWrapper& df, const CITestType& ci_test_type);
