#pragma once
#include <cstdint>
#include <cstdio>
#include <vector>

#include "base/contingency_table.h"
#include "citest/citest_type.h"

/*
 * @brief Perform a conditional independence test between two variables x and y
 * given a conditioning set Z.
 * Time complexity: O(|X| * |Y| * ∏_{Z_i ∈ Z} |Z_i|) where |X|, |Y| are the
 * cardinalities of x and y and each |Z_i| is the cardinality of a variable in
 * Z.
 *
 * @param x The first variable index.
 * @param y The second variable index.
 * @param sepset_candidate The conditioning set Z.
 * @param ct The contingency counts for the variables.
 * @param ci_test_type The type of CI test to perform.
 * @return true if x and y are conditionally independent given Z, false
 * otherwise.
 */
bool citest(size_t x, size_t y, const std::vector<size_t>& sepset_candidate,
            const ContingencyTable& ct, const CITestType& ci_test_type);

// std::vector<std::vector<int>> state_count(const
// std::vector<std::vector<uint8_t>> &data,
//                                 std::vector<int> &children, std::vector<int>
//                                 &parents, std::vector<int> &n_states, int
//                                 &parallel);

// float natori_independent_score(const std::vector<std::vector<uint8_t>> &data,
// int &node_x,
//                                int &node_y, std::vector<int> &parents,
//                                std::vector<int> &n_states, int &parallel);

// float natori_dependent_score(const std::vector<std::vector<uint8_t>> &data,
// int &node_x,
//                              int &node_y, std::vector<int> &parents,
//                              std::vector<int> &n_states, int &parallel);

// bool ci_test_by_Bayesfactor(const std::vector<std::vector<uint8_t>> &data,
// int &node_x,
//                             int &node_y, std::vector<int> &Z,
//                             std::vector<int> &n_states, int &parallel);