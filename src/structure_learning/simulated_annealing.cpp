#include "structure_learning/simulated_annealing.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

#include "score/score_type.h"
#include "score/parent_set_evaluator.h"
#include "utils/logging.h"
#include "utils/timeout_guard.h"

PDAG simulated_annealing(const DataframeWrapper& df,
                         const ScoreType& score_type,
                         size_t max_parents,
                         size_t max_iters,
                         double init_temp,
                         double cooling_rate,
                         uint64_t seed,
                         size_t num_chains,
                         double timeout_sec) {
  double MINUS_INF = -std::numeric_limits<double>::infinity();
  TimeoutGuard tg(timeout_sec);
  ParentSetEvaluator pse(df, score_type, max_parents, tg);

  if (num_chains == 0) num_chains = omp_get_max_threads();
  INFO("[SA] running " << num_chains << " chains in parallel ...");

  const size_t n = df.num_vars;
  std::vector<double> best_score_of_chain(num_chains, MINUS_INF);
  std::vector<std::vector<std::vector<size_t>>> best_parents_of_chain(num_chains);

#pragma omp parallel for schedule(dynamic)
  for (size_t chain_id = 0; chain_id < num_chains; ++chain_id) {
    std::mt19937_64 rng(seed + chain_id);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    
    double global_score = 0.0;
    std::vector<double> local_scores(n, 0.0);
    std::vector<std::vector<size_t>> parents(n);

    // Initialize a DAG without edges, i.e., no cycles
    for (size_t child = 0; child < n; ++child) {
      const auto& parent_set_cands = pse.p_pars[child];
      for (const auto& [ls, pars] : parent_set_cands) {
        if (pars.size() == 0) {
          global_score += ls;
          local_scores[child] = ls;
          parents[child] = pars;
          break;
        }
      }
    }

    double best_global_score = global_score;
    auto best_local_scores = local_scores;
    auto best_parents = parents;
    double T = init_temp;
    for (size_t it = 0; it < max_iters; ++it) { // for each chain

      T *= cooling_rate;
      if (T < 1e-12) T = 1e-12;

      // Propose a new state by changing the parent set of a randomly chosen node
      size_t child = rng() % n;
      const auto& parent_set_cands = pse.p_pars[child];
      size_t new_pars_idx = rng() % parent_set_cands.size();
      const auto& [new_score, new_pars] = parent_set_cands[new_pars_idx];

      // check acyclicity
      bool creates_cycle = false;
      std::vector<bool> visited(n, false);
      std::vector<size_t> stack;
      for (const auto& p : new_pars) {
        visited[p] = true;
        stack.push_back(p);
      }
      while (!stack.empty() && !creates_cycle) {
        size_t node = stack.back();
        stack.pop_back();
        for (const auto& parent : parents[node]) {
          if (parent == child) {
            creates_cycle = true;
            break;
          }
          if (!visited[parent]) {
            visited[parent] = true;
            stack.push_back(parent);
          }
        }
      }
      if (creates_cycle) continue;

      // delta score
      double delta = new_score - local_scores[child];
      bool accept = (delta >= 0.0) || (unif(rng) < std::exp(delta / T));
      if (accept) {
        global_score += delta;
        local_scores[child] = new_score;
        parents[child] = new_pars;
        if (global_score > best_global_score) {
          best_global_score = global_score;
          best_local_scores = local_scores;
          best_parents = parents;
        }
      }

      if (it % 10000 == 0) {
        if (tg.is_timeout() || tg.expired()) break;  // check timeout
      }
    }

    best_score_of_chain[chain_id] = best_global_score;
    best_parents_of_chain[chain_id] = best_parents;
  }

  size_t best_idx = 0;
  for (size_t c = 1; c < best_score_of_chain.size(); ++c) {
    if (best_score_of_chain[c] > best_score_of_chain[best_idx]) best_idx = c;
  }

  PDAG best_dag(df.num_vars);
  for (size_t child = 0; child < df.num_vars; ++child) {
    for (const auto& parent : best_parents_of_chain[best_idx][child]) {
      best_dag.add_edge(parent, child);
    }
  }
  return best_dag;
}
