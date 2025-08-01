#include "structure_learning/exhaustive_search.h"

#include <iostream>

#include "base/contingency_table.h"
#include "score/local_score.h"
#include "utils/comb2vec.h"
#include "utils/next_combination.h"

#define varset_t uint64_t  // number of variables <= 64

PDAG exhaustive_search(const DataframeWrapper& df,
                       const ScoreType& score_type,
                       size_t max_parents,
                       bool is_deterministic) {
  if (max_parents < 0 || max_parents > df.num_of_vars - 1)
    throw std::invalid_argument("max_parents must be in [0, num_of_vars-1]");
  size_t max_varset_size = max_parents + 1;
  size_t n = df.num_of_vars;
  if (n < 0 || n > 64)
    throw std::invalid_argument("num_of_vars must be in [0, 64]");

  std::vector<varset_t> bit_masks(n);
  for (size_t i = 0; i < n; ++i) {
    bit_masks[i] = (varset_t(1) << i);
  }

  size_t num_of_combs = 1 << n;

  // variable set score table in size 2^n
  std::vector<double> best_gs_tbl(num_of_combs);
  std::vector<uint8_t> best_ch_tbl(num_of_combs);

  // best local score and its parents table in size n * 2^n
  std::vector<std::vector<std::pair<double, varset_t>>> best_ls_tbl(n);
  for (size_t i = 0; i < n; i++) {
    best_ls_tbl[i].resize(num_of_combs);  // todo devide by 2
  }

  double empty_score = 0;
  varset_t empty_varset_int = 0;
  best_gs_tbl[empty_varset_int] = empty_score;

  // explore all possible variable subsets in size from 0 to n
  for (size_t subset_size = 0; subset_size <= n; subset_size++) {
    bool calc_ls = subset_size <= max_varset_size;

    // explore all possible variable combinations in the current subset size
    varset_t varset_int = (varset_t(1) << subset_size) - 1;
    do {
      std::vector<size_t> varset_vec = comb2vec(varset_int);

      // TODO
      // ContingencyTable<true> ct;
      // if (calc_ls) ct = ContingencyTable<true>(varset_vec, df);
      ContingencyTable<true> ct(varset_vec, df);

      double best_gs = 0;
      int best_ch = -1;
      for (size_t child_idx = 0; child_idx < varset_vec.size(); child_idx++) {
        auto child_var = varset_vec[child_idx];
        auto parents_varset_int = varset_int - bit_masks[child_var];

        double ls = -std::numeric_limits<double>::infinity();
        if (calc_ls) {
          std::vector<size_t> parent_set;
          for (auto v : varset_vec) {
            if (v != child_var) parent_set.push_back(v);
          }
          ls = calculate_local_score<double>(
              child_var, parent_set, ct, score_type);
        }

        double best_ls = ls;
        varset_t best_ps = parents_varset_int;
        for (const auto& parent_var : varset_vec) {
          if (parent_var == child_var) continue;
          auto parent_subset_int = parents_varset_int - bit_masks[parent_var];
          auto [best_ls_of_subset, best_ps_of_subset] =
              best_ls_tbl[child_var][parent_subset_int];
          if (best_ls_of_subset >= best_ls) {
            best_ls = best_ls_of_subset;
            best_ps = best_ps_of_subset;
          }
        }
        best_ls_tbl[child_var][parents_varset_int] = {best_ls, best_ps};

        double gs = best_gs_tbl[parents_varset_int] + best_ls;
        // since best_gs is always negative, 0 means uninitialized
        if (best_gs >= 0 || gs > best_gs) {
          best_gs = gs;
          best_ch = child_var;
        }

        // /* debug print begin */
        // std::cout << "child: " << df.col_idx2str[child_var]
        //           << ", parents: [";
        // for (const auto& parent : varset_vec) {
        //   if (parent == child_var) continue;
        //   std::cout << df.col_idx2str[parent] << ", ";
        // }
        // std::cout << "]"
        //           << ", ls: " << ls
        //           << ", best_ls: " << best_ls
        //           << ", best_ps: [";
        // for (size_t j = 0; j < n; j++) {
        //   if (best_ps & bit_masks[j]) {
        //     std::cout << df.col_idx2str[j] << ", ";
        //   }
        // }
        // std::cout << "]" << std::endl;
        // /* debug print end */
      }
      best_gs_tbl[varset_int] = best_gs;
      best_ch_tbl[varset_int] = best_ch;

      // /* debug print begin */
      // std::cout << "varset: [";
      // for (const auto& var : varset_vec) {
      //   std::cout << df.col_idx2str[var] << ", ";
      // }
      // std::cout << "]"
      //           << ", best_gs: " << best_gs
      //           << ", best_ch: " << df.col_idx2str[best_ch]
      //           << std::endl;
      // /* debug print end */

    } while (next_combination(varset_int, n));
  }

  // reconstruct the best graph
  PDAG best_graph(n);
  varset_t varset_int = (varset_t(1) << n) - 1;
  for (size_t i = 0; i < n; i++) {
    // auto best_gs = best_gs_tbl[varset_int];
    auto best_ch = best_ch_tbl[varset_int];

    auto parents_varset_int = varset_int - bit_masks[best_ch];
    auto [best_ls, best_ps] = best_ls_tbl[best_ch][parents_varset_int];

    for (size_t j = 0; j < n; j++) {
      if (best_ps & bit_masks[j]) {
        best_graph.add_edge(j, best_ch);
      }
    }
    varset_int = parents_varset_int;
  }

  // /* debug print begin */
  // std::cout << "best gs: " << best_gs_tbl[(varset_t(1) << n) - 1]
  //           << std::endl;
  // for (size_t i = 0; i < n; i++) {
  //   std::cout << df.col_idx2str[i] << " <- {";
  //   for (size_t j = 0; j < n; j++) {
  //     if (best_graph.has_edge(j, i)) {
  //       std::cout << df.col_idx2str[j] << ", ";
  //     }
  //   }
  //   std::cout << "}" << std::endl;
  // }
  // /* debug print end */

  return best_graph;
}