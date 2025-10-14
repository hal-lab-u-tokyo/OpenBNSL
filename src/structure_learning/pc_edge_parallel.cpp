#include "structure_learning/pc_edge_parallel.h"

#include "base/contingency_table.h"
#include "base/dataframe_wrapper.h"
#include "citest/citest.h"
#include "graph/pdag_with_adjmat.h"
#include "utils/gen_comb.h"
#include "utils/logging.h"
#include "utils/timeout_guard.h"

struct Update {
  size_t x, y;
  std::vector<size_t> Z;
};
namespace pc_edge_parallel_ns {
void build_skeleton(PDAGwithAdjMat& g,
                    Sepset& sepset,
                    const DataframeWrapper& df,
                    const CITestType& test,
                    size_t max_cond_vars,
                    TimeoutGuard& tg,
                    std::vector<size_t>& num_citests_summary,
                    size_t& exit_order) {
  const size_t n = g.num_vars;
  for (size_t k = 0; k <= max_cond_vars; ++k) {
    if (tg.expired()) throw std::runtime_error("Timeout");
    auto start = std::chrono::high_resolution_clock::now();

    // Precompute neighbors
    std::vector<std::vector<size_t>> neighs(n);
    for (size_t x = 0; x < n; ++x) {
      neighs[x] = g.undirected_neighbors(x);
    }

    /* Exit condition */
    bool has_sepset_candidates = false;
    for (size_t x = 0; x < n; ++x) {
      if (neighs[x].size() >= k + 1) {
        has_sepset_candidates = true;
        break;
      }
    }
    if (!has_sepset_candidates) {
      exit_order = k;
      break;
    }

    std::vector<std::pair<size_t, size_t>> pairs;
    for (size_t x = 0; x < n; ++x) {
      for (size_t y = x + 1; y < n; ++y) {
        if (g.has_undirected_edge(x, y)) {
          pairs.emplace_back(x, y);
        }
      }
    }

    size_t num_citests = 0;
    std::vector<Update> updates;
#pragma omp parallel default(none) \
    shared(n, k, df, test, pairs, updates, neighs, num_citests, tg)
    {
      size_t local_num_citests = 0;
      std::vector<Update> local_updates;
#pragma omp for schedule(dynamic)
      for (size_t i = 0; i < pairs.size(); ++i) {
        if (tg.is_timeout() || tg.expired()) continue;  // check timeout
        auto [x, y] = pairs[i];
        // iterate on smaller one first
        if (neighs[x].size() > neighs[y].size()) std::swap(x, y);
        for (auto [u, v] : std::array{std::pair{x, y}, std::pair{y, x}}) {
          const auto& neigh_of_u = neighs[u];
          if (neigh_of_u.size() < k + 1) continue;
          std::vector<size_t> neigh_of_u_without_v;
          neigh_of_u_without_v.reserve(neigh_of_u.size() - 1);
          for (auto w : neigh_of_u)
            if (w != v) neigh_of_u_without_v.push_back(w);
          std::sort(neigh_of_u_without_v.begin(), neigh_of_u_without_v.end());

          for (const auto& Z : gen_combs(neigh_of_u_without_v, k)) {
            local_num_citests++;
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
      }  // end omp for
#pragma omp critical
      {
        num_citests += local_num_citests;
        updates.insert(
            updates.end(), local_updates.begin(), local_updates.end());
      }
    }  // end omp parallel

    for (const auto& update : updates) {
      g.remove_edge(update.x, update.y);
      sepset[update.x][update.y].insert(update.Z.begin(), update.Z.end());
      sepset[update.y][update.x] = sepset[update.x][update.y];
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    INFO("[PC] order=" << k << ", edges removed=" << updates.size() << "/"
                       << num_citests << ", time=" << elapsed.count() << "s");
    num_citests_summary[k] += num_citests;
  }
}
}  // namespace pc_edge_parallel_ns

PDAG pc_edge_parallel(const DataframeWrapper& df,
                      const CITestType& test,
                      size_t max_cond_vars,
                      double timeout_sec) {
  // Only for logging
  size_t exit_order = 0;
  std::vector<size_t> num_citests_summary(max_cond_vars + 1, 0);

  TimeoutGuard tg(timeout_sec);
  const size_t n = df.num_vars;
  PDAGwithAdjMat g(n);
  g.set_as_complete();
  Sepset sepset(n, std::vector<std::unordered_set<size_t>>(n));
  pc_edge_parallel_ns::build_skeleton(
      g, sepset, df, test, max_cond_vars, tg, num_citests_summary, exit_order);
  g.orient_colliders(sepset);
  g.apply_meeks_rules();

  // Only for logging
  auto vec_to_str = [](const std::vector<size_t>& vec, size_t max_len = 0) {
    std::string s = "[";
    for (size_t i = 0; i < vec.size() && i < max_len; ++i) {
      s += std::to_string(i) + ":" + std::to_string(vec[i]) + ", ";
    }
    s += "]";
    return s;
  };
  INFO("[PC] total CI tests:" << vec_to_str(num_citests_summary, exit_order));

  return g.to_pdag();
}
