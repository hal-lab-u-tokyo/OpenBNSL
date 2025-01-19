#pragma once
#include <optional>
#include <vector>

#include "score/score_type.h"

/**
 * @brief Calculate the local score
 * @tparam T Return type (float or double)
 * @param child_var Child variable (target for score calculation)
 * @param varset_vec Variable set vector, including child_var
 * @param freq_tbl Frequency table for calculations
 * @param score_type Score type (default: BDeu)
 * @param ess Equivalent Sample Size (optional, required for BDeu)
 * @return The calculated local score as the specified type T
 * @throw std::invalid_argument If score_type is unsupported or required
 * parameters are missing
 */
template <typename T>
T calculate_local_score(int child_var,
                        const std::vector<int>& varset,  // including child_var
                        const std::vector<size_t>& varset_size,
                        const std::vector<int>& freq_tbl,
                        const ScoreType& score_type);
