#include "structure_learning/pc.h"

#include <algorithm>
#include <array>
#include <vector>

#include "base/contingency_table.h"
#include "base/dataframe_wrapper.h"
#include "citest/citest.h"
#include "graph/pdag_with_adjmat.h"
#include "utils/gen_comb.h"

void build_skeleton(PDAGwithAdjMat& g,
                    const DataframeWrapper& df,
                    const CITestType& test,
                    size_t max_cond_vars,
                    Sepset& sepset) {
  const size_t n = g.num_global_vars;

  for (size_t k = 0; k <= max_cond_vars; ++k) {
    const PDAGwithAdjMat snapshot = g;
    struct Update {
      size_t x, y;
      std::vector<size_t> Z;
    };
    std::vector<Update> updates;

#pragma omp parallel default(none) shared(n, df, test, snapshot, k, updates)
    {
      std::vector<Update> local_updates;
#pragma omp for schedule(dynamic)
      for (size_t x = 0; x < n; ++x) {
        for (size_t y = x + 1; y < n; ++y) {
          if (!snapshot.has_undirected_edge(x, y)) continue;

          for (auto [u, v] : std::array{std::pair{x, y}, std::pair{y, x}}) {
            auto neigh = snapshot.undirected_neighbors_without(u, v);
            if (neigh.size() < k) continue;

            // B: Create contingency table for {u,v} U neigh(u) once
            // std::vector<size_t> vars = neigh;
            // vars.push_back(u);
            // vars.push_back(v);
            // std::sort(vars.begin(), vars.end());
            // ContingencyTable ct(vars, df);

            for (const auto& Z : gen_combs(neigh, k)) {
              // A: Create contingency table for {u,v} U Z for
              std::vector<size_t> vars = Z;
              vars.push_back(u);
              vars.push_back(v);
              std::sort(vars.begin(), vars.end());
              ContingencyTable ct(vars, df);

              if (citest(u, v, Z, ct, test)) {
                Update update{x, y, std::vector<size_t>(Z.begin(), Z.end())};
                local_updates.push_back(update);
                goto NEXT_PAIR;
              }
            }
          }
        NEXT_PAIR:;
        }
      }
#pragma omp critical
      updates.insert(updates.end(), local_updates.begin(), local_updates.end());
    }  // end omp parallel

    for (const auto& update : updates) {
      g.remove_undirected_edge(update.x, update.y);
      sepset[update.x][update.y].insert(update.Z.begin(), update.Z.end());
      sepset[update.y][update.x] = sepset[update.x][update.y];
    }
  }
}

PDAG pc(const DataframeWrapper& df,
        const CITestType& test,
        size_t max_cond_vars) {
  const size_t n = df.num_of_vars;
  PDAGwithAdjMat g(n);
  g.set_as_complete();
  Sepset sepset(n, std::vector<std::unordered_set<size_t>>(n));

  build_skeleton(g, df, test, max_cond_vars, sepset);
  g.orient_colliders(sepset);
  g.apply_meeks_rules();
  return g.to_pdag();
}
