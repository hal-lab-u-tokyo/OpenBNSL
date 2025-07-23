#include "structure_learning/simulated_annealing.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

#include "base/contingency_table.h"
#include "graph/pdag_with_adjlist.h"
#include "score/local_score.h"

template <typename GraphT>
static void run_single_chain(const DataframeWrapper& df,
                             const ScoreType& score_type,
                             int max_parents,
                             size_t max_iters,
                             double init_temp,
                             double cooling_rate,
                             uint64_t seed,
                             double& best_score_out,
                             GraphT& best_graph_out) {
  const size_t n = df.num_of_vars;
  GraphT g(n);
  std::vector<double> ls(n, 0.0);
  for (size_t v = 0; v < n; ++v) {
    ls[v] = calculate_local_score<double>(
        v, {}, buildContingencyTable({v}, df), score_type);
  }

  double curr_score = std::accumulate(ls.begin(), ls.end(), 0.0);
  double best_score = curr_score;
  auto best_g = g;

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  std::uniform_int_distribution<size_t> pick_node(0, n - 1);

  double T = init_temp;
  for (size_t it = 0; it < max_iters; ++it, T *= cooling_rate) {
    size_t child = pick_node(rng);
    const auto& pa = g.parents[child];
    enum { ADD, REMOVE } op;

    if (pa.empty())
      op = ADD;
    else if (pa.size() == static_cast<size_t>(max_parents))
      op = REMOVE;
    else
      op = (uni(rng) < 0.5 ? ADD : REMOVE);

    std::vector<size_t> new_pa(pa.begin(), pa.end());

    if (op == ADD) {
      std::vector<size_t> cand(n);
      std::iota(cand.begin(), cand.end(), 0);
      cand.erase(std::remove(cand.begin(), cand.end(), child), cand.end());
      std::shuffle(cand.begin(), cand.end(), rng);

      bool done = false;
      for (auto p : cand) {
        if (pa.count(p)) continue;
        if (g.has_path(p, child)) continue;
        new_pa.push_back(p);
        std::sort(new_pa.begin(), new_pa.end());
        done = true;
        break;
      }
      if (!done) continue;
    } else {
      std::uniform_int_distribution<size_t> pick_par(0, new_pa.size() - 1);
      new_pa.erase(new_pa.begin() + pick_par(rng));
    }

    std::vector<size_t> vars = new_pa;
    vars.push_back(child);
    std::sort(vars.begin(), vars.end());
    double new_ls = calculate_local_score<double>(
        child, new_pa, buildContingencyTable(vars, df), score_type);
    double delta = new_ls - ls[child];

    if (delta >= 0.0 || std::exp(delta / T) > uni(rng)) {
      g.set_parents(
          child, typename GraphT::ParentSetType(new_pa.begin(), new_pa.end()));
      ls[child] = new_ls;
      curr_score += delta;
      if (curr_score > best_score) {
        best_score = curr_score;
        best_g = g;
      }
    }
  }

  best_score_out = best_score;
  best_graph_out = std::move(best_g);
}

PDAG simulated_annealing(const DataframeWrapper& df,
                         const ScoreType& score_type,
                         int max_parents,
                         size_t max_iters,
                         double init_temp,
                         double cooling_rate,
                         bool is_deterministic,
                         uint64_t seed,
                         int num_chains) {
  if (max_parents < 0 || max_parents >= static_cast<int>(df.num_of_vars)) {
    throw std::invalid_argument("max_parents out of range");
  }

  if (num_chains <= 0) {
    num_chains = omp_get_max_threads();
  }

  std::vector<double> scores(num_chains, -1e100);
  PDAG result(df.num_of_vars);

  if (is_deterministic) {
    using GraphT = PDAGwithAdjList<true>;
    std::vector<GraphT> graphs(num_chains, GraphT(df.num_of_vars));
#pragma omp parallel for
    for (int c = 0; c < num_chains; ++c) {
      run_single_chain<GraphT>(df,
                               score_type,
                               max_parents,
                               max_iters,
                               init_temp,
                               cooling_rate,
                               seed + c * 1234567ULL,
                               scores[c],
                               graphs[c]);
    }
    int best_idx = std::distance(
        scores.begin(), std::max_element(scores.begin(), scores.end()));
    result = graphs[best_idx].to_pdag();
  } else {
    using GraphT = PDAGwithAdjList<false>;
    std::vector<GraphT> graphs(num_chains, GraphT(df.num_of_vars));
#pragma omp parallel for
    for (int c = 0; c < num_chains; ++c) {
      run_single_chain<GraphT>(df,
                               score_type,
                               max_parents,
                               max_iters,
                               init_temp,
                               cooling_rate,
                               seed + c * 1234567ULL,
                               scores[c],
                               graphs[c]);
    }
    int best_idx = std::distance(
        scores.begin(), std::max_element(scores.begin(), scores.end()));
    result = graphs[best_idx].to_pdag();
  }

  return result;
}
