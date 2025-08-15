#pragma once

#include <algorithm>
#include <boost/math/distributions/chi_squared.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "base/contingency_table.h"
#include "citest/citest_type.h"
#include "citest/graph_dsep.h"
#include "graph/pdag.h"

inline double pchisq(double x, std::size_t dof) {
  if (dof == 0) return (x == 0.0 ? 1.0 : 0.0);
  const boost::math::chi_squared dist(dof);
  return 1.0 - boost::math::cdf(dist, x);
}

/**
 * @ingroup citest
 * @brief Conditional–independence test for discrete variables.
 * @details
 * The contingency table @p ct is stored sparsely, i.e. cells with observed
 * count 0 do not appear.  pgmpy/Scipy, however, still accounts for such cells
 * if the expected frequency is non‑zero.  We therefore generate those cells
 * on‑the‑fly: for each Z‑slice we iterate over the Cartesian product of the
 * actually observed X–states and Y–states and fill in obs = 0 where necessary.
 * @tparam Deterministic Whether the contingency table is deterministic or not
 * @param x The first variable
 * @param y The second variable
 * @param sepset_candidate The candidate separator set
 * @param ct The contingency table
 * @param ci_test_type The CI test type
 * @return True if the conditional independence test passes (p_value >= alpha),
 *         false otherwise.
 */
template <bool Deterministic>
bool citest(std::size_t x,
            std::size_t y,
            const std::vector<std::size_t>& sepset_candidate,
            const ContingencyTable<Deterministic>& ct,
            const CITestType& ci_test_type) {
  if (x == y) throw std::invalid_argument("x and y must differ");
  if (x >= y) std::swap(x, y);

  double alpha;
  if (is_type<ChiSquare>(ci_test_type)) {
    alpha = static_cast<double>(get_type<ChiSquare>(ci_test_type).level);
  } else if (is_type<GSquare>(ci_test_type)) {
    alpha = static_cast<double>(get_type<GSquare>(ci_test_type).level);
  } else if (is_type<OracleGraph>(ci_test_type)) {
    const PDAG& g = get_type<OracleGraph>(ci_test_type).graph;
    return is_d_separated(g, x, y, sepset_candidate);
  } else {
    throw std::invalid_argument("Unsupported CI test type");
  }

  const auto it_x = std::find(ct.var_ids.begin(), ct.var_ids.end(), x);
  const auto it_y = std::find(ct.var_ids.begin(), ct.var_ids.end(), y);
  if (it_x == ct.var_ids.end() || it_y == ct.var_ids.end())
    throw std::invalid_argument("x or y not contained in contingency table");

  const std::size_t idx_x = static_cast<std::size_t>(it_x - ct.var_ids.begin());
  const std::size_t idx_y = static_cast<std::size_t>(it_y - ct.var_ids.begin());

  const std::size_t n_x = ct.cardinalities[idx_x];
  const std::size_t n_y = ct.cardinalities[idx_y];

  // Generate slices for each observed z state
  std::unordered_map<std::size_t, std::vector<std::size_t>> slices;
  for (const auto& [key, freq] : ct.counts) {
    const std::size_t x_state = ct.state_of(key, idx_x);
    const std::size_t y_state = ct.state_of(key, idx_y);
    std::size_t slice_key = key - x_state * ct.radix_weight(idx_x) -
                            y_state * ct.radix_weight(idx_y);

    auto& slice = slices[slice_key];
    if (slice.empty()) slice.assign(n_x * n_y, 0);
    slice[x_state * n_y + y_state] += freq;
  }

  double stat = 0;
  std::size_t dof = 0;
  for (const auto& [_, slice] : slices) {
    std::vector<double> row_sum(n_x, 0), col_sum(n_y, 0);
    double n_total = 0;
    for (std::size_t i = 0; i < n_x; ++i) {
      for (std::size_t j = 0; j < n_y; ++j) {
        const double obs = static_cast<double>(slice[i * n_y + j]);
        row_sum[i] += obs;
        col_sum[j] += obs;
        n_total += obs;
      }
    }

    if (n_total == 0) continue;  // skip empty slices
    std::vector<std::size_t> rows_obs, cols_obs;
    for (std::size_t i = 0; i < n_x; ++i) {
      if (row_sum[i] > 0) rows_obs.push_back(i);
    }
    for (std::size_t j = 0; j < n_y; ++j) {
      if (col_sum[j] > 0) cols_obs.push_back(j);
    }

    const std::size_t _dof = (rows_obs.size() - 1) * (cols_obs.size() - 1);
    const bool apply_yates = (_dof == 1);
    dof += _dof;

    for (auto i : rows_obs) {
      for (auto j : cols_obs) {
        double obs = slice[i * n_y + j];
        double exp = row_sum[i] * col_sum[j] / n_total;

        // Apply Yates' correction (reference:
        // https://github.com/scipy/scipy/blob/v1.16.1/scipy/stats/contingency.py)
        if (apply_yates) {
          double diff = exp - obs;
          double direction = (diff > 0) ? 1.0 : (diff < 0) ? -1.0 : 0.0;
          double magnitude = std::min(0.5, std::fabs(diff));
          obs += magnitude * direction;
        }

        if (is_type<ChiSquare>(ci_test_type)) {
          double diff = obs - exp;
          stat += diff * diff / exp;
        } else if (is_type<GSquare>(ci_test_type)) {
          if (obs > 0) stat += 2.0 * obs * std::log(obs / exp);
        } else {
          throw std::invalid_argument("Unsupported CI test type");
        }
      }
    }
  }

  const double p_value = pchisq(stat, dof);
  // std::cout << "[obnsl] stat: " << stat
  //           << ", p_value: " << p_value
  //           << ", dof: " << dof
  //           << ", level: " << alpha
  //           << ", result: " << (p_value >= alpha ? "True" : "False") <<
  //           std::endl;

  return p_value >= alpha;
}
