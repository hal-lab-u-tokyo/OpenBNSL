#include "structure_learning/PC.h"

#include <cstddef>
#include <vector>

#include "base/contingency_table.h"
#include "base/dataframe_wrapper.h"
#include "citest/citest.h"
#include "utils/next_combination.h"

PDAG build_skeleton(const DataframeWrapper& df, const CITestType& citest_type) {
  size_t n = df.num_of_vars;
  PDAG g(n);
  g.complete_graph();  // Start with a complete graph

  size_t l = 0;  // the order of the conditioning set
  while (l <= n - 2) {
    for (size_t x = 0; x < n; ++x) {
      std::vector<size_t> x_neighbors = g.neighbors(x);
      // TODO: In case not enough neighbors to condition on

      std::vector<size_t> varset(x_neighbors);
      varset.push_back(x);
      std::sort(varset.begin(), varset.end());
      ContingencyTable ct = buildContingencyTable(varset, df);

      for (size_t y : x_neighbors) {
        // Get the neighbors of x excluding y
        std::vector<size_t> adj_without_y;
        for (size_t neighbor : x_neighbors) {
          if (neighbor != y) {
            adj_without_y.push_back(neighbor);
          }
        }

        // Not enough neighbors to condition on
        if (adj_without_y.size() < l) continue;

        // Generate all combinations of size l from adj_without_y
        std::vector<bool> comb_mask(adj_without_y.size(), false);
        std::fill(comb_mask.begin(), comb_mask.begin() + l, true);
        do {
          std::vector<size_t> sepset_candidate;
          for (size_t i = 0; i < comb_mask.size(); ++i) {
            if (comb_mask[i]) {
              sepset_candidate.push_back(adj_without_y[i]);
            }
          }
          if (citest(x, y, sepset_candidate, ct, citest_type)) {
            g.remove_edge(x, y);
            g.remove_edge(y, x);
            break;  // No need to check further sepset candidates
          }
        } while (std::prev_permutation(comb_mask.begin(), comb_mask.end()));
      }
    }
    l++;
  }

  return g;
}

PDAG PC2(const DataframeWrapper& df, const CITestType& citest_type) {
  PDAG g = build_skeleton(df, citest_type);

  return g;
}
