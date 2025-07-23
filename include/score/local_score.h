#pragma once
#include <vector>

#include "base/contingency_table.h"
#include "score/score_type.h"

/**
 * @brief Calculate the local score
 * @tparam T Return type (float or double)
 * @param child_var Child variable
 * @param parent_set Parent variable set
 * @param ct Contingency table
 * @param score_type Score type (default: BDeu)
 * @return The calculated local score as the specified precision type
 */
template <typename T>
T calculate_local_score(size_t child_var,
                        const std::vector<size_t>& parent_set,
                        const ContingencyTable& ct,
                        const ScoreType& score_type = BDeu{1.0});
