#pragma once
#include <cstdint>
#include <cstdio>
#include <vector>

std::vector<uint8_t> debug_int(std::vector<uint8_t> a);

// std::vector<std::vector<int>> state_count(const std::vector<std::vector<uint8_t>> &data,
//                                 std::vector<int> &children, std::vector<int> &parents,
//                                 std::vector<int> &n_states, int &parallel);

// float natori_independent_score(const std::vector<std::vector<uint8_t>> &data, int
// &node_x,
//                                int &node_y, std::vector<int> &parents,
//                                std::vector<int> &n_states, int &parallel);

// float natori_dependent_score(const std::vector<std::vector<uint8_t>> &data, int
// &node_x,
//                              int &node_y, std::vector<int> &parents,
//                              std::vector<int> &n_states, int &parallel);

// bool ci_test_by_Bayesfactor(const std::vector<std::vector<uint8_t>> &data, int &node_x,
//                             int &node_y, std::vector<int> &Z, std::vector<int>
//                             &n_states, int &parallel);