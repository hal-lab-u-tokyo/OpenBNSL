#include "citest/citest.h"

#include <algorithm>
#include <boost/math/distributions/chi_squared.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "citest/graph_dsep.h"
#include "graph/pdag.h"
#include "utils/enumerate.h"

double pchisq(double x, std::size_t dof) {
  if (dof == 0) return (x == 0.0 ? 1.0 : 0.0);
  const boost::math::chi_squared dist(dof);
  return 1.0 - boost::math::cdf(dist, x);
}

bool citest(std::size_t x,
            std::size_t y,
            const std::vector<std::size_t>& sepset_candidate,
            const ContingencyTable& ct,
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

  std::size_t idx_x = ct.var_ids.size();
  std::size_t idx_y = ct.var_ids.size();
  std::size_t itr = 0;
  std::vector<std::size_t> idx_extras;
  for (auto [idx, var_id] : enumerate(ct.var_ids)) {
    if (var_id == x)
      idx_x = idx;
    else if (var_id == y)
      idx_y = idx;
    else if (itr < sepset_candidate.size() && var_id == sepset_candidate[itr])
      ++itr;
    else
      idx_extras.push_back(idx);
  }
  if (idx_x == ct.var_ids.size() || idx_y == ct.var_ids.size())
    throw std::invalid_argument("x or y not contained in contingency table");
  if (itr != sepset_candidate.size())
    throw std::invalid_argument(
        "sepset candidate not contained in contingency table");

  const std::size_t n_x = ct.cardinalities[idx_x];
  const std::size_t n_y = ct.cardinalities[idx_y];

  // Generate slices for each observed z state
  std::unordered_map<std::size_t, std::unordered_map<std::size_t, std::size_t>>
      slices_sparse;
  for (const auto& [key, freq] : ct.counts) {
    const std::size_t x_state = ct.state_of(key, idx_x);
    const std::size_t y_state = ct.state_of(key, idx_y);

    std::size_t slice_key = key;
    slice_key = ct.strip(slice_key, idx_x);
    slice_key = ct.strip(slice_key, idx_y);
    for (auto idx : idx_extras) slice_key = ct.strip(slice_key, idx);

    std::size_t xy_key = x_state * n_y + y_state;
    slices_sparse[slice_key][xy_key] += freq;
  }

  double stat = 0;
  std::size_t dof = 0;
  for (const auto& [sepset_key, slice] : slices_sparse) {
    std::unordered_map<std::size_t, std::size_t> row_sum, col_sum;
    double n_total = 0;
    for (const auto& [xy_key, freq] : slice) {
      const std::size_t x_state = xy_key / n_y;
      const std::size_t y_state = xy_key % n_y;
      row_sum[x_state] += freq;
      col_sum[y_state] += freq;
      n_total += freq;
    }

    const std::size_t _dof = (row_sum.size() - 1) * (col_sum.size() - 1);
    const bool apply_yates = (_dof == 1);
    dof += _dof;

    for (const auto& [xy_key, obs_raw] : slice) {
      const std::size_t x_state = xy_key / n_y;
      const std::size_t y_state = xy_key % n_y;
      double obs = static_cast<double>(obs_raw);
      double exp = static_cast<double>(row_sum[x_state]) *
                   static_cast<double>(col_sum[y_state]) /
                   static_cast<double>(n_total);

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

  const double p_value = pchisq(stat, dof);
  // std::cout << "[obnsl] stat: " << stat
  //           << ", p_value: " << p_value
  //           << ", dof: " << dof
  //           << ", level: " << alpha
  //           << ", result: " << (p_value >= alpha ? "True" : "False") <<
  //           std::endl;

  return p_value >= alpha;
}