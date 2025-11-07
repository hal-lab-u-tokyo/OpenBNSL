
#include "score/parent_set_evaluator.h"

#include <algorithm>
#include <chrono>
#include <numeric>
#include <stdexcept>

#include "score/local_score.h"
#include "utils/logging.h"
#include "utils/binom.h"
#include "utils/combvec.h"
#include "utils/vector_utils.h"

ParentSetEvaluator::ParentSetEvaluator(const DataframeWrapper& df_,
                 const ScoreType& score_type_,
                 size_t max_parents_,
                 TimeoutGuard& tg)
                 : df(df_), score_type(score_type_), max_parents(max_parents_) {
  auto start = std::chrono::high_resolution_clock::now();
                
  size_t n = df.num_vars;
  if (max_parents >= n) max_parents = n - 1;
  nCk_tbl = utils::build_binom_table(n, max_parents);
  // TODO: more sophisticated check
  // if (true) throw std::runtime_error("Problem is too large");
  
  ls_tbl.resize(n);
  for (size_t x = 0; x < n; ++x) {
    ls_tbl[x].resize(max_parents + 1);
    for (size_t k = 0; k <= max_parents; ++k) {
      ls_tbl[x][k].resize(nCk_tbl[n - 1][k], 0);
    }
  }
  p_pars.resize(n);

  std::vector<size_t> all_vars(n);
  std::iota(all_vars.begin(), all_vars.end(), 0);
#pragma omp parallel for schedule(dynamic)
  for (size_t x = 0; x < n; ++x) {
    if (tg.expired()) throw std::runtime_error("Timeout");
    std::vector<size_t> vars_without_x = utils::copy_without(all_vars, x);
    auto ct = ContingencyTable({x}, df);
    ls_tbl[x][0][0] = calculate_local_score(x, {}, ct, score_type);
    p_pars[x].emplace_back(ls_tbl[x][0][0], std::vector<size_t>{});
    for (size_t k = 1; k <= max_parents; ++k) {
      if (tg.expired()) throw std::runtime_error("Timeout");
      for (const auto& Z : utils::gen_combvecs(vars_without_x, k)) {
        std::vector<size_t> vars = Z;
        vars.push_back(x);
        std::sort(vars.begin(), vars.end());
        ContingencyTable ct(vars, df);
        double score = calculate_local_score(x, Z, ct, score_type);
        size_t rank = utils::combvec_rank(utils::reindex_excluding(x, Z), nCk_tbl);
        ls_tbl[x][k][rank] = score;
        bool prune = false;
        for (auto y : Z) {
          auto subset = utils::copy_without(Z, y);
          size_t subset_rank = utils::combvec_rank(utils::reindex_excluding(x, subset), nCk_tbl);
          if (ls_tbl[x][k - 1][subset_rank] >= score) {
            prune = true;
            break;
          }
        }
        if (!prune) p_pars[x].emplace_back(score, Z);
      }
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  auto summary = [&](const std::vector<std::vector<std::pair<double, std::vector<size_t>>>>& all) {
    std::string s = "[";
    for (size_t x = 0; x < all.size(); ++x) {
      s += std::to_string(all[x].size());
      if (x + 1 < all.size()) s += ", ";
    }
    s += "]";
    return s;
  };
  INFO("[ParentSetEvaluator] ready. #cands per var = " << summary(p_pars)
       << ", time = " << elapsed.count() << "s");
}

double& ParentSetEvaluator::at(size_t x, size_t k, const std::vector<size_t>& S) {
  auto idx = utils::reindex_excluding(x, S);
  size_t rank = utils::combvec_rank(idx, nCk_tbl);
  return ls_tbl.at(x).at(k).at(rank);
}

const double& ParentSetEvaluator::at(size_t x,
                          size_t k,
                          const std::vector<size_t>& S) const {
  auto idx = utils::reindex_excluding(x, S);
  size_t rank = utils::combvec_rank(idx, nCk_tbl);
  return ls_tbl.at(x).at(k).at(rank);
}

