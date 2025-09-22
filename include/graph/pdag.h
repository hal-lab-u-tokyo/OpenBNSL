#pragma once
#include <cstddef>
#include <set>
#include <vector>

#include "base/contingency_table.h"
#include "base/dataframe_wrapper.h"
#include "score/local_score.h"
#include "score/score_type.h"

/**
 * @ingroup graph
 * @struct PDAG
 * @brief Interface for converting internal graph structures to a unified PDAG.
 * @details
 * Structure learning algorithms may use different graph representations.
 * This interface ensures that all such graphs can be converted into a common
 * sparse PDAG form for downstream use (e.g., evaluation, comparison, or Python
 * interop).
 */
struct PDAG {
  /* Data members */
  std::size_t num_vars;
  std::vector<std::set<size_t>> parents;  // parents[v] = {u | u -> v}

  /* Lifecycle */
  /**
   * @brief Construct a new PDAG with a specified number of variables.
   * @param num_vars The number of variables in the PDAG.
   */
  PDAG(std::size_t num_vars) : num_vars(num_vars), parents(num_vars) {}

  bool has_edge(std::size_t from, std::size_t to) const {
    return parents[to].find(from) != parents[to].end();
  }

  void add_edge(std::size_t from, std::size_t to) { parents[to].insert(from); }

  void remove_edge(std::size_t from, std::size_t to) {
    parents[to].erase(from);
  }

  double score(const DataframeWrapper& df, const ScoreType& score_type) const {
    double res = 0.0;
    for (std::size_t v = 0; v < num_vars; ++v) {
      std::vector<size_t> parents_vec;
      for (auto pa : parents[v]) parents_vec.push_back(pa);
      std::vector<size_t> vars = parents_vec;
      vars.push_back(v);
      std::sort(vars.begin(), vars.end());
      ContingencyTable<false> ct(vars, df);
      res +=
          calculate_local_score<double, false>(v, parents_vec, ct, score_type);
    }
    return res;
  }
};
