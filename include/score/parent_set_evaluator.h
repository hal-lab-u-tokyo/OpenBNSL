#pragma once
#include <cstddef>
#include <vector>

#include "base/dataframe_wrapper.h"
#include "score/score_type.h"
#include "utils/timeout_guard.h"

/**
 * @brief Local score cache for local scores
 */
struct ParentSetEvaluator {
  const DataframeWrapper& df;
  const ScoreType& score_type;
  size_t max_parents;

  std::vector<std::vector<size_t>> nCk_tbl;

  // ls_tbl[x][k][rank]:
  // Local score ls(x, S) for the k-sized parent set S \subset {0,...,n-1}
  // \ {x}. The parent set S is identified by its combination rank `rank`.
  std::vector<std::vector<std::vector<double>>> ls_tbl;

  // p_pars[x][i]:
  // The i-th best parent set candidate for variable x,
  // represented as a pair (score, parent_set).
  std::vector<std::vector<std::pair<double, std::vector<size_t>>>> p_pars;

  explicit ParentSetEvaluator(const DataframeWrapper& df,
                              const ScoreType& score_type,
                              size_t max_parents,
                              TimeoutGuard& tg);
  double& at(size_t x, size_t k, const std::vector<size_t>& S);
  const double& at(size_t x, size_t k, const std::vector<size_t>& S) const;
};
