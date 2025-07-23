#pragma once
#include "base/dataframe_wrapper.h"
#include "graph/pdag.h"
#include "score/score_type.h"

/**
 * @brief Score-based multi-chain simulated annealing
 * @param df Dataframe containing the data
 * @param score_type Type of score to use (e.g., BIC, AIC)
 * @param max_parents Maximum number of parents for each node
 * @param max_iters Maximum number of iterations for the annealing process
 * @param init_temp Initial temperature for the annealing process
 * @param cooling_rate Cooling rate for the temperature
 * @param seed Random seed for reproducibility
 * @param num_chains Number of parallel chains to run (0 for auto-detect)
 * @param is_deterministic Whether to use deterministic or non-deterministic
 * @return Best DAG found during the annealing process
 */
PDAG simulated_annealing(const DataframeWrapper& df,
                         const ScoreType& score_type,
                         int max_parents,
                         size_t max_iters = 10000,
                         double init_temp = 1.0,
                         double cooling_rate = 0.9995,
                         bool is_deterministic = false,
                         uint64_t seed = 0,
                         int num_chains = 0);
