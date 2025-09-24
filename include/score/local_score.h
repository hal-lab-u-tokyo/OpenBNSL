#pragma once

#include <cstddef>
#include <vector>

#include "base/contingency_table.h"
#include "score/score_type.h"

/**
 * @ingroup score
 * @brief Calculate the local score for a given child variable.
 * @param child_var The child variable for which to calculate the score.
 * @param parent_set The set of parent variables.
 * @param ct The contingency table.
 * @param score_type The type of score to calculate.
 * @return The calculated local score.
 */
double calculate_local_score(size_t child_var,
                             const std::vector<size_t>& /*parent_set*/,
                             const ContingencyTable& ct,
                             const ScoreType& score_type = BDeu{1.0});