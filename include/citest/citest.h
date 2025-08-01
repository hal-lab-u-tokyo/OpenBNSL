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

template <typename StatScalar>
inline StatScalar pchisq(StatScalar x, std::size_t dof) {
  if (dof == 0) return (x == 0.0 ? 1.0 : 0.0);
  const boost::math::chi_squared dist(dof);
  return 1.0 - boost::math::cdf(dist, x);
}

/**
 * @brief Conditional–independence test for discrete variables.
 *
 * The contingency table @p ct is stored sparsely, i.e. cells with observed
 * count 0 do not appear.  pgmpy/Scipy, however, still accounts for such cells
 * if the expected frequency is non‑zero.  We therefore generate those cells
 * on‑the‑fly: for each Z‑slice we iterate over the Cartesian product of the
 * actually observed X–states and Y–states and fill in obs = 0 where necessary.
 *
 * @return True if the conditional independence test passes (p_value >= alpha),
 *         false otherwise.
 */
template <typename StatScalar, bool Deterministic>
bool citest(std::size_t                       x,
            std::size_t                       y,
            const std::vector<std::size_t>&   /*sepset_candidate*/,
            const ContingencyTable<Deterministic>& ct,
            const CITestType&                 ci_test_type)
{
  if (x >= y) std::swap(x, y);

  StatScalar alpha;
  if (is_type<ChiSquare>(ci_test_type))
    alpha = static_cast<StatScalar>(get_type<ChiSquare>(ci_test_type).level);
  else if (is_type<GSquare>(ci_test_type))
    alpha = static_cast<StatScalar>(get_type<GSquare>(ci_test_type).level);
  else
    throw std::invalid_argument("Unsupported CI test type");

  const auto it_x = std::find(ct.var_ids.begin(), ct.var_ids.end(), x);
  const auto it_y = std::find(ct.var_ids.begin(), ct.var_ids.end(), y);
  if (it_x == ct.var_ids.end() || it_y == ct.var_ids.end())
    throw std::invalid_argument("x or y not contained in contingency table");

  const std::size_t idx_x = static_cast<std::size_t>(it_x - ct.var_ids.begin());
  const std::size_t idx_y = static_cast<std::size_t>(it_y - ct.var_ids.begin());

  const std::size_t rw_x  = ct.radix_weight(idx_x);
  const std::size_t rw_y  = ct.radix_weight(idx_y);

  CountsMap<Deterministic> n_xz, n_yz, n_z;
  
  // unique X and Y values for each Z slice  
  std::unordered_map<std::size_t, std::unordered_set<std::size_t>> x_vals_of_z;
  std::unordered_map<std::size_t, std::unordered_set<std::size_t>> y_vals_of_z;

  for (const auto& [key, cnt] : ct.counts) {
    const std::size_t key_xz = ct.strip(key, idx_y);
    const std::size_t key_yz = ct.strip(key, idx_x);
    const std::size_t key_z  = ct.strip(key_xz, idx_x);

    n_xz[key_xz] += cnt;
    n_yz[key_yz] += cnt;
    n_z [key_z ] += cnt;

    x_vals_of_z[key_z].insert(ct.state_of(key, idx_x));
    y_vals_of_z[key_z].insert(ct.state_of(key, idx_y));
  }

  StatScalar stat = 0.0;
  for (const auto& [key_z, pz] : n_z) {

    std::cerr << "[Z] key=" << key_z << " N=" << pz
            << "  |X|=" << xs_of_z[key_z].size()
            << "  |Y|=" << ys_of_z[key_z].size() << '\n';

    const auto& xs = x_vals_of_z[key_z];
    const auto& ys = y_vals_of_z[key_z];

    for (std::size_t xv : xs) {
      const std::size_t key_xz = key_z + xv * rw_x;
      const StatScalar  pxz    = static_cast<StatScalar>(n_xz[key_xz]);

      for (std::size_t yv : ys) {
        const std::size_t key_yz = key_z + yv * rw_y;
        const StatScalar  pyz    = static_cast<StatScalar>(n_yz[key_yz]);

        const StatScalar  exp = pxz * pyz / static_cast<StatScalar>(pz);

        const std::size_t full_key = key_z + xv * rw_x + yv * rw_y;
        const StatScalar  obs = ct.contains(full_key)
                                ? static_cast<StatScalar>(ct.counts.at(full_key))
                                : static_cast<StatScalar>(0);

                                
        if (is_type<ChiSquare>(ci_test_type)) {
          const StatScalar diff = obs - exp;
          stat += (diff * diff) / exp;
        } else {
          if (obs > 0)
            stat += static_cast<StatScalar>(2.0) * obs * std::log(obs / exp);
        }
      }
    }
  }

  std::size_t dof = 0;
  for (const auto& [key_z, pz] : n_z) {
    if (pz == 0) continue;
    const std::size_t r = x_vals_of_z[key_z].size();
    const std::size_t c = y_vals_of_z[key_z].size();
    if (r >= 2 && c >= 2) dof += (r - 1) * (c - 1);
  }

  const StatScalar p_value = pchisq(stat, dof);
  std::cout << "[obnsl] stat: " << stat 
            << ", p_value: " << p_value
            << ", dof: " << dof 
            << ", level: " << alpha
            << ", result: " << (p_value >= alpha ? "True" : "False") << std::endl;

  return p_value >= alpha;
}
