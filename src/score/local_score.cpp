#include "score/local_score.h"

#include <cmath>
#include <stdexcept>

#include "score/score_type.h"

template <typename T>
T calculate_local_score(int child_var, const std::vector<int>& varset,
                        const std::vector<size_t>& varset_size,
                        const std::vector<int>& freq_tbl,
                        const ScoreType& score_type) {
  size_t outer_size = 1, child_size = 1, inner_size = 1;
  for (const auto& var : varset) {
    if (var < child_var) {
      outer_size *= varset_size[var];
    } else if (var == child_var) {
      child_size = varset_size[var];
    } else {  // var > child_var
      inner_size *= varset_size[var];
    }
  }
  if (outer_size * child_size * inner_size != freq_tbl.size()) {
    throw std::invalid_argument("Invalid frequency table size.");
  }

  T a_ijk;
  // using ScoreType = std::variant<BDeu, K2>;
  // switch (score_type) {
  //   case ScoreType::BDeu:
  //     a_ijk = ess.value_or(1.0) / (outer_size * child_size * inner_size);
  //     break;
  //   case ScoreType::K2:
  //     a_ijk = 1.0;
  //     break;
  //   default:
  //     throw std::invalid_argument("Unsupported score type.");
  // }
  if (is_type<BDeu>(score_type)) {
    const auto& bdeu = get_type<BDeu>(score_type);
    a_ijk = (T)bdeu.ess / (outer_size * child_size * inner_size);
  } else if (is_type<K2>(score_type)) {
    a_ijk = 1.0;
  } else {
    throw std::invalid_argument("Unsupported score type.");
  }
  T a_ij = a_ijk * child_size;
  T lgamma_a_ij = std::lgamma(a_ij);
  T lgamma_a_ijk = std::lgamma(a_ijk);

  T ls = 0;
  size_t inner_idx = 0;
  for (size_t i = 0; i < outer_size; i++) {
    std::vector<int> N_ij(inner_size, 0);
    for (size_t j = 0; j < child_size; j++) {
      for (size_t k = 0; k < inner_size; k++) {
        // sequentially access the frequency table
        // to improve cache locality and hit rate
        auto N_ijk = freq_tbl[inner_idx++];
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
    int child_var, const std::vector<int>& varset,
    const std::vector<size_t>& varset_size, const std::vector<int>& freq_tbl,
    const ScoreType& score_type);

template double calculate_local_score<double>(
    int child_var, const std::vector<int>& varset,
    const std::vector<size_t>& varset_size, const std::vector<int>& freq_tbl,
    const ScoreType& score_type);
