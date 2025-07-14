#include "score/local_score.h"

#include <cmath>
#include <stdexcept>

#include "score/score_type.h"

template <typename T>
T calculate_local_score(size_t child_var, const std::vector<size_t>& parent_set,
                        const ContingencyTable& ct,
                        const ScoreType& score_type) {
  if ( is_type<BIC>(score_type) || is_type<AIC>(score_type) ) {
    throw std::invalid_argument("Unsupported score type.");
  }
  
  size_t outer_size = 1, child_size = 1, inner_size = 1;
  for (size_t i = 0; i < ct.var_ids.size(); i++) {
    if (ct.var_ids[i] < child_var) {
      outer_size *= ct.cardinalities[i];
    } else if (ct.var_ids[i] == child_var) {
      child_size = ct.cardinalities[i];
    } else {  // ct.var_ids[i] > child_var
      inner_size *= ct.cardinalities[i];
    }
  }
  size_t total_size = outer_size * child_size * inner_size;
  if (total_size != ct.counts.size()) {
    throw std::invalid_argument("Invalid frequency table size.");
  }

  T a_ijk;
  if (is_type<BDeu>(score_type)) {
    const auto& bdeu = get_type<BDeu>(score_type);
    a_ijk = (T)bdeu.ess / total_size;
  } else if (is_type<K2>(score_type)) {
    a_ijk = 1.0;
  } else {
    throw std::invalid_argument("Unsupported score type.");
  }
  T a_ij = a_ijk * child_size;
  T lgamma_a_ij = std::lgamma(a_ij);
  T lgamma_a_ijk = std::lgamma(a_ijk);

  T ls = 0;
  size_t addr = 0;
  for (size_t i = 0; i < outer_size; i++) {
    std::vector<int> N_ij(inner_size, 0);
    for (size_t j = 0; j < child_size; j++) {
      for (size_t k = 0; k < inner_size; k++) {
        auto N_ijk = ct.counts[addr++];
        ls += std::lgamma(N_ijk + a_ijk) - lgamma_a_ijk;
        N_ij[k] += N_ijk;
      }
    }
    for (size_t k = 0; k < inner_size; k++) {
      ls -= std::lgamma(N_ij[k] + a_ij) - lgamma_a_ij;
    }
  }

  return ls;
}

template float calculate_local_score<float>(
    size_t child_var, const std::vector<size_t>& parent_set,
    const ContingencyTable& ct, const ScoreType& score_type);

template double calculate_local_score<double>(
    size_t child_var, const std::vector<size_t>& parent_set,
    const ContingencyTable& ct, const ScoreType& score_type);
