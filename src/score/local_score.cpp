#include "score/local_score.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

double calculate_local_score(size_t child_var,
                             const std::vector<size_t>& /*parent_set*/,
                             const ContingencyTable& ct,
                             const ScoreType& score_type) {
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

  double a_ijk;
  if (is_type<BDeu>(score_type)) {
    const auto& bdeu = get_type<BDeu>(score_type);
    a_ijk = (double)bdeu.ess / total_size;
  } else {
    throw std::invalid_argument("Unsupported score type.");
  }
  double a_ij = a_ijk * child_card;
  double lgamma_a_ijk = std::lgamma(a_ijk);
  double lgamma_a_ij = std::lgamma(a_ij);

  double ls = 0;
  std::unordered_map<size_t, size_t> marged_counts;
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