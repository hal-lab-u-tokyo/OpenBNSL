#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
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
                                  const ScoreType& score_type = BDeu{1.0}) {
  auto itr = std::find(ct.var_ids.begin(), ct.var_ids.end(), child_var);
  if (itr == ct.var_ids.end()) {
    throw std::invalid_argument("child_var not present in contingency table");
  }
  const size_t child_idx = std::distance(ct.var_ids.begin(), itr);
  const size_t child_card = ct.cardinalities[child_idx];
  size_t total_size = 1;
  for (auto card : ct.cardinalities) {
    total_size *= card;
  }

  ScoreScalar a_ijk;
  if (is_type<BDeu>(score_type)) {
    const auto& bdeu = get_type<BDeu>(score_type);
    a_ijk = (ScoreScalar)bdeu.ess / total_size;
  } else {
    throw std::invalid_argument("Unsupported score type.");
  }
  ScoreScalar a_ij = a_ijk * child_card;
  ScoreScalar lgamma_a_ijk = std::lgamma(a_ijk);
  ScoreScalar lgamma_a_ij = std::lgamma(a_ij);

  ScoreScalar ls = 0;
  CountsMap<Deterministic> marged_counts;
  for (const auto& [key, N_ijk] : ct.counts) {
    if (N_ijk == 0) continue;
    ls += std::lgamma(N_ijk + a_ijk) - lgamma_a_ijk;
    const size_t parent_key = ct.strip(key, child_idx);
    marged_counts[parent_key] += N_ijk;
  }
  for (const auto& [_, N_ij] : marged_counts) {
    if (N_ij == 0) continue;
    ls -= std::lgamma(N_ij + a_ij) - lgamma_a_ij;
  }
  return ls;
}