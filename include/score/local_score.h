#pragma once

#include <cstddef>
#include <vector>

#include "base/contingency_table.h"
#include "score/score_type.h"

/**
 * @ingroup score
 * @brief Calculate the local score for a given child variable.
 * @tparam ScoreScalar The type used for the score calculation.
 * @tparam Deterministic Whether the calculation is deterministic.
 * @param child_var The child variable for which to calculate the score.
 * @param parent_set The set of parent variables.
 * @param ct The contingency table.
 * @param score_type The type of score to calculate.
 * @return The calculated local score.
 */
template <typename ScoreScalar, bool Deterministic>
ScoreScalar calculate_local_score(size_t child_var,
                                  const std::vector<size_t>& /*parent_set*/,
                                  const ContingencyTable<Deterministic>& ct,
                                  const ScoreType& score_type = BDeu{1.0});