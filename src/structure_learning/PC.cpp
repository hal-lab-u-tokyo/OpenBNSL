#include "structure_learning/pc.h"

#include <algorithm>
#include <array>
#include <unordered_set>
#include <vector>

#include "base/contingency_table.h"
#include "base/dataframe_wrapper.h"
#include "citest/citest.h"
#include "graph/pdag_with_adjmat.h"
#include "utils/gen_comb.h"

using Sepset = std::vector<std::vector<std::unordered_set<size_t>>>;

static PDAGwithAdjMat build_skeleton(const DataframeWrapper& df,
                                     const CITestType& test,
                                     size_t max_cond_vars,
                                     bool stable,
                                     Sepset& sepset) {
  PDAGwithAdjMat g(df.num_of_vars);
  for (size_t k = 0; k <= max_cond_vars; ++k) {

    // snapshot the current graph if stable version
    PDAGwithAdjMat snapshot = g;
    PDAGwithAdjMat* g_ref = stable ? &snapshot : &g;

    for (size_t x = 0; x < g.num_vars; ++x) {
      for (size_t y = x + 1; y < g.num_vars; ++y) {
        // Skip if there is no edge between x and y
        if (!g.has_undirected_edge(x, y)) continue;
        for (auto [u, v] : std::array{std::pair{x, y}, std::pair{y, x}}) {
          auto neigh = g_ref->undirected_neighbors_without(u, v);
          for (auto& Z : gen_combs(neigh, k)) {
            std::vector<size_t> vars = Z;
            vars.push_back(u);
            vars.push_back(v);
            std::sort(vars.begin(), vars.end());
            ContingencyTable<true> ct(vars, df);
            if (citest<true>(u, v, Z, ct, test)) {
              g.remove_undirected_edge(x, y);
              sepset[x][y].insert(Z.begin(), Z.end());
              sepset[y][x] = sepset[x][y];
              goto NEXT_PAIR;
            }
          }
        }
      NEXT_PAIR:;
      }
    }
  }
  return g;
}

void orient_colliders(PDAGwithAdjMat& g, const Sepset& sepset) {
  const size_t n = g.num_vars;
  for (size_t y = 0; y < n; ++y) {
    bool changed = true;
    while (changed) {
      changed = false;
      auto neigh = g.undirected_neighbors(y);
      for (size_t i = 0; i + 1 < neigh.size(); ++i) {
        for (size_t j = i + 1; j < neigh.size(); ++j) {
          size_t x = neigh[i], z = neigh[j];
          if (g.is_adjacent(x, z)) continue;
          if (sepset[x][z].count(y) == 0) {
            g.orient_edge(x, y);
            g.orient_edge(z, y);
            changed = true;
          }
        }
      }
    }
  }
}

PDAG pc(const DataframeWrapper& df,
        const CITestType& test,
        size_t max_cond_vars,
        bool stable) {
  const size_t n = df.num_of_vars;
  Sepset sepset(n, std::vector<std::unordered_set<size_t>>(n));
  PDAGwithAdjMat g = build_skeleton(df, test, max_cond_vars, stable, sepset);
  orient_colliders(g, sepset);
  // g.apply_meeks_rules(false);
  return g.to_pdag();
}
