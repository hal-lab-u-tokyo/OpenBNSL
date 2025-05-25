#include "citest/citest.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <cmath>
#include <boost/math/distributions/chi_squared.hpp>

double pchisq(double x, size_t dof) {
    if (dof == 0) return 1.0;
    boost::math::chi_squared dist(dof);
    return 1.0 - boost::math::cdf(dist, x);
}

bool citest(
  size_t x, size_t y, 
  const std::vector<size_t>& sepset_candidate,
  const ContingencyTable& ct, const CITestType& ci_test_type) {

  double alpha;
  if (is_type<ChiSquare>(ci_test_type)) {
      alpha = get_type<ChiSquare>(ci_test_type).level;
  }
  else if (is_type<GSquare>(ci_test_type)) {
      alpha = get_type<GSquare>(ci_test_type).level;
  }
  else {
      throw std::invalid_argument("Unsupported CI test type");
  }
  
    // Ensure x < y
  if (x >= y) std::swap(x, y);
  // Ensure sepset_candidate is sorted
  for (size_t i = 1; i < sepset_candidate.size(); ++i) {
    if (sepset_candidate[i] < sepset_candidate[i - 1]) {
      throw std::invalid_argument("sepset_candidate must be sorted");
    }
  }

  // Split sepset_candidate into three parts: z0, z1, z2
  // Ensure _z0 < x < _z1 < y < _z2,
  // , where forall _z0 \in z0, _z1 \in z1, _z2 \in z2
  // , and z0 \sqcup z1 \sqcup z2 = sepset_candidate
  size_t n_x = 1, n_y = 1, n_z0 = 1, n_z1 = 1, n_z2 = 1;
  for (size_t i = 0; i < ct.var_ids.size(); ++i) {
      if      (ct.var_ids[i] <  x) n_z0 *= ct.cardinalities[i];
      else if (ct.var_ids[i] == x) n_x   = ct.cardinalities[i];
      else if (ct.var_ids[i] <  y) n_z1 *= ct.cardinalities[i];
      else if (ct.var_ids[i] == y) n_y   = ct.cardinalities[i];
      else                         n_z2 *= ct.cardinalities[i];
  }
  if (n_z0 * n_x * n_z1 * n_y * n_z2 != ct.counts.size())
      throw std::invalid_argument("Contingency table size mismatch");

  size_t addr = 0;  // Ensure the sequential access
  size_t total_size_x = n_z0 * n_z1 * n_z2 * n_x;
  size_t total_size_y = n_z0 * n_z1 * n_z2 * n_y;
  size_t total_size_xy = n_z0 * n_z1 * n_z2;
  std::vector<size_t> marginalized_x(total_size_x, 0);
  std::vector<size_t> marginalized_y(total_size_y, 0);
  std::vector<size_t> marginalized_xy(total_size_xy, 0);
  auto index_x = [=](size_t idx_z0, size_t idx_z1, size_t idx_z2,
                     size_t idx_x) {
    return (idx_z0 * n_z1 + idx_z1) * n_z2 * n_x + idx_z2 * n_x + idx_x;
  };
  auto index_y = [=](size_t idx_z0, size_t idx_z1, size_t idx_z2,
                     size_t idx_y) {
    return (idx_z0 * n_z1 + idx_z1) * n_z2 * n_y + idx_z2 * n_y + idx_y;
  };
  auto index_xy = [=](size_t idx_z0, size_t idx_z1, size_t idx_z2) {
    return (idx_z0 * n_z1 + idx_z1) * n_z2 + idx_z2;
  };
  double stat = 0.0;

  addr = 0;
  for (size_t idx_z0 = 0; idx_z0 < n_z0; ++idx_z0) {
    for (size_t idx_x = 0; idx_x < n_x; ++idx_x) {
      for (size_t idx_z1 = 0; idx_z1 < n_z1; ++idx_z1) {
        for (size_t idx_y = 0; idx_y < n_y; ++idx_y) {
          for (size_t idx_z2 = 0; idx_z2 < n_z2; ++idx_z2) {
            size_t observed_count = ct.counts[addr++];
            marginalized_x[index_x(idx_z0, idx_z1, idx_z2, idx_x)] +=
                observed_count;
            marginalized_y[index_y(idx_z0, idx_z1, idx_z2, idx_y)] +=
                observed_count;
            marginalized_xy[index_xy(idx_z0, idx_z1, idx_z2)] += observed_count;
          }
        }
      }
    }
  }

  addr = 0;
  for (size_t idx_z0 = 0; idx_z0 < n_z0; ++idx_z0) {
    for (size_t idx_x = 0; idx_x < n_x; ++idx_x) {
      for (size_t idx_z1 = 0; idx_z1 < n_z1; ++idx_z1) {
        for (size_t idx_y = 0; idx_y < n_y; ++idx_y) {
          for (size_t idx_z2 = 0; idx_z2 < n_z2; ++idx_z2) {
            double observed_count = static_cast<double>(ct.counts[addr++]);
            double denom = marginalized_xy[index_xy(idx_z0, idx_z1, idx_z2)];
            if (denom == 0) continue;  // Avoid division by zero
            double expected_count =
                static_cast<double>(
                    marginalized_x[index_x(idx_z0, idx_z1, idx_z2, idx_x)] *
                    marginalized_y[index_y(idx_z0, idx_z1, idx_z2, idx_y)]) / denom;
            if (expected_count == 0) continue;
            
            if (is_type<ChiSquare>(ci_test_type)) {
              double diff = observed_count - expected_count;
              stat += (diff * diff) / expected_count;
            } else if (is_type<GSquare>(ci_test_type)) {
              if (observed_count == 0.0) continue;  // Avoid log(0)
              stat += 2.0 * observed_count * std::log(observed_count / expected_count);
            } else {
              throw std::invalid_argument("Unsupported CI test type");
            }
          }
        }
      }
    }
  }

  size_t dof = (n_x - 1) * (n_y - 1) * n_z0 * n_z1 * n_z2;
  double p_value = pchisq(stat, dof);
  return p_value > alpha; // true if independent, false if dependent
}