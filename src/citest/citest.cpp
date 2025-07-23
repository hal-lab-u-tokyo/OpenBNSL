#include "citest/citest.h"

#include <algorithm>
#include <boost/math/distributions/chi_squared.hpp>
#include <cmath>
#include <functional>
#include <iostream>

double pchisq(double x, std::size_t dof) {
  if (dof == 0) return 1.0;
  boost::math::chi_squared dist(dof);
  return 1.0 - boost::math::cdf(dist, x);
}

bool citest(std::size_t x,
            std::size_t y,
            const std::vector<std::size_t>& sepset_candidate,
            const ContingencyTable& ct,
            const CITestType& ci_test_type) {
  double alpha;
  if (is_type<ChiSquare>(ci_test_type)) {
    alpha = get_type<ChiSquare>(ci_test_type).level;
  } else if (is_type<GSquare>(ci_test_type)) {
    alpha = get_type<GSquare>(ci_test_type).level;
  } else {
    throw std::invalid_argument("Unsupported CI test type");
  }

  // Ensure x < y
  if (x >= y) std::swap(x, y);
  // Ensure sepset_candidate is sorted
  for (std::size_t i = 1; i < sepset_candidate.size(); ++i) {
    if (sepset_candidate[i] < sepset_candidate[i - 1]) {
      throw std::invalid_argument("sepset_candidate must be sorted");
    }
  }

  // Split sepset_candidate into three parts: z0, z1, z2
  // Ensure _z0 < x < _z1 < y < _z2,
  // , where forall _z0 \in z0, _z1 \in z1, _z2 \in z2
  // , and z0 \sqcup z1 \sqcup z2 = sepset_candidate
  std::size_t n_x = 1, n_y = 1, n_z0 = 1, n_z1 = 1, n_z2 = 1;
  for (std::size_t i = 0; i < ct.var_ids.size(); ++i) {
    if (ct.var_ids[i] < x)
      n_z0 *= ct.cardinalities[i];
    else if (ct.var_ids[i] == x)
      n_x = ct.cardinalities[i];
    else if (ct.var_ids[i] < y)
      n_z1 *= ct.cardinalities[i];
    else if (ct.var_ids[i] == y)
      n_y = ct.cardinalities[i];
    else
      n_z2 *= ct.cardinalities[i];
  }
  if (n_z0 * n_x * n_z1 * n_y * n_z2 != ct.counts.size())
    throw std::invalid_argument("Contingency table size mismatch");

  // Create indices for accessing counts in the contingency table
  auto index_xz = [=](std::size_t i_z0,
                      std::size_t i_x,
                      std::size_t i_z1,
                      std::size_t i_z2) {
    return ((i_z0 * n_x + i_x) * n_z1 + i_z1) * n_z2 + i_z2;
  };
  auto index_yz = [=](std::size_t i_z0,
                      std::size_t i_z1,
                      std::size_t i_y,
                      std::size_t i_z2) {
    return ((i_z0 * n_z1 + i_z1) * n_y + i_y) * n_z2 + i_z2;
  };
  auto index_z = [=](std::size_t i_z0, std::size_t i_z1, std::size_t i_z2) {
    return (i_z0 * n_z1 + i_z1) * n_z2 + i_z2;
  };

  std::size_t addr = 0;  // Ensure the sequential access
  std::vector<std::size_t> p_xz(n_z0 * n_x * n_z1 * n_z2, 0);
  std::vector<std::size_t> p_yz(n_z0 * n_z1 * n_y * n_z2, 0);
  std::vector<std::size_t> p_z(n_z0 * n_z1 * n_z2, 0);
  addr = 0;
  for (std::size_t i_z0 = 0; i_z0 < n_z0; ++i_z0) {
    for (std::size_t i_x = 0; i_x < n_x; ++i_x) {
      for (std::size_t i_z1 = 0; i_z1 < n_z1; ++i_z1) {
        for (std::size_t i_y = 0; i_y < n_y; ++i_y) {
          for (std::size_t i_z2 = 0; i_z2 < n_z2; ++i_z2) {
            const std::size_t obs_cnt = ct.counts[addr++];
            p_xz[index_xz(i_z0, i_x, i_z1, i_z2)] += obs_cnt;
            p_yz[index_yz(i_z0, i_z1, i_y, i_z2)] += obs_cnt;
            p_z[index_z(i_z0, i_z1, i_z2)] += obs_cnt;
          }
        }
      }
    }
  }

  std::vector<bool> valid_z(n_z0 * n_z1 * n_z2, true);
  for (std::size_t i_z0 = 0; i_z0 < n_z0; ++i_z0) {
    for (std::size_t i_x = 0; i_x < n_x; ++i_x) {
      for (std::size_t i_z1 = 0; i_z1 < n_z1; ++i_z1) {
        for (std::size_t i_y = 0; i_y < n_y; ++i_y) {
          for (std::size_t i_z2 = 0; i_z2 < n_z2; ++i_z2) {
            std::size_t _p_xz = p_xz[index_xz(i_z0, i_x, i_z1, i_z2)];
            std::size_t _p_yz = p_yz[index_yz(i_z0, i_z1, i_y, i_z2)];
            if (_p_xz == 0 || _p_yz == 0) {
              valid_z[index_z(i_z0, i_z1, i_z2)] = false;
            }
          }
        }
      }
    }
  }

  std::size_t vld_cnt = 0;
  for (auto v : valid_z) {
    if (v) ++vld_cnt;
  }

  double stat = 0.0;
  addr = 0;
  for (std::size_t i_z0 = 0; i_z0 < n_z0; ++i_z0) {
    for (std::size_t i_x = 0; i_x < n_x; ++i_x) {
      for (std::size_t i_z1 = 0; i_z1 < n_z1; ++i_z1) {
        for (std::size_t i_y = 0; i_y < n_y; ++i_y) {
          for (std::size_t i_z2 = 0; i_z2 < n_z2; ++i_z2) {
            if (valid_z[index_z(i_z0, i_z1, i_z2)]) {
              double obs_cnt = static_cast<double>(ct.counts[addr]);
              double _p_xz =
                  static_cast<double>(p_xz[index_xz(i_z0, i_x, i_z1, i_z2)]);
              double _p_yz =
                  static_cast<double>(p_yz[index_yz(i_z0, i_z1, i_y, i_z2)]);
              double _p_z = static_cast<double>(p_z[index_z(i_z0, i_z1, i_z2)]);
              double exp_cnt = _p_xz * _p_yz / _p_z;

              if (is_type<ChiSquare>(ci_test_type)) {
                double diff = obs_cnt - exp_cnt;
                stat += (diff * diff) / exp_cnt;
              } else if (is_type<GSquare>(ci_test_type)) {
                if (obs_cnt == 0.0 || exp_cnt < 1e-8) continue;
                stat += 2.0 * obs_cnt * std::log(obs_cnt / exp_cnt);
              } else {
                throw std::invalid_argument("Unsupported CI test type");
              }
            }
            ++addr;
          }
        }
      }
    }
  }
  std::size_t dof = (n_x - 1) * (n_y - 1) * vld_cnt;

  double p_value = pchisq(stat, dof);
  std::cout << "[obnel] stat: " << stat << ", p_value: " << p_value
            << ", dof: " << dof << ", alpha: " << alpha
            << ", result: " << (p_value >= alpha) << std::endl;

  return p_value >= alpha;  // true if independent, false if dependent
}