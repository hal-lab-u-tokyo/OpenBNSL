#include "structure_learning/PC.h"

#include <algorithm>
#include <queue>
#include <set>
#include <unordered_set>

#include "base/contingency_table.h"
#include "base/dataframe_wrapper.h"
#include "citest/citest.h"

void enumerate_combinations(
  const std::vector<size_t> &items, 
  size_t k,
  size_t idx, std::vector<size_t> &current,
  std::vector<std::vector<size_t>> &result
) {
  if (current.size() == k) {
    result.push_back(current);
    return;
  }
  for (size_t i = idx; i < items.size(); ++i) {
    current.push_back(items[i]);
    enumerate_combinations(items, k, i + 1, current, result);
    current.pop_back();
  }
}

std::vector<std::vector<size_t>> combinations(
  const std::vector<size_t> &items,
  size_t k) {
  std::vector<std::vector<size_t>> res;
  if (k == 0) {
    res.push_back({});
    return res;
  }
  std::vector<size_t> current;
  enumerate_combinations(items, k, 0, current, res);
  return res;
}

static bool all_nodes_degree_less_than(const PDAG &g, size_t limit) {
    const size_t n = g.num_vars;
    for (size_t v = 0; v < n; ++v) {
        if (g.undirected_neighbors(v).size() >= limit) return false;
    }
    return true;
}

PDAG build_skeleton(
  const DataframeWrapper &df, 
  const CITestType &citest_type,
  size_t max_cond_vars,
  std::vector<std::vector<std::vector<size_t>>> &sepset
) {
  const size_t n = df.num_of_vars;
  PDAG pdag(n);
  pdag.complete_graph();
  
  for (size_t l = 0; l <= max_cond_vars; ++l) {
    for (size_t x = 0; x < n; ++x) {
      for (size_t y = x + 1; y < n; ++y) {
        if (!pdag.has_edge(x, y)) continue; // already removed
        std::vector<size_t> neigh_x = pdag.undirected_neighbors(x);
        neigh_x.erase(std::remove(neigh_x.begin(), neigh_x.end(), y), neigh_x.end());
        for (auto &sepset_candidates: combinations(neigh_x, l)) {
          std::vector<size_t> var_ids(sepset_candidates);
          var_ids.push_back(x);
          var_ids.push_back(y);
          std::sort(var_ids.begin(), var_ids.end());
          ContingencyTable ct = buildContingencyTable(var_ids, df);
          if (citest(x, y, sepset_candidates, ct, citest_type)) {
            pdag.remove_edge(x, y);
            pdag.remove_edge(y, x);
            sepset[x][y] = sepset[y][x] = sepset_candidates;
            break;  // stop checking other sets for this pair
          }
        } 

        if (!pdag.has_edge(x, y)) continue; // already removed
        std::vector<size_t> neigh_y = pdag.undirected_neighbors(y);
        neigh_y.erase(std::remove(neigh_y.begin(), neigh_y.end(), x), neigh_y.end());
        for (auto &sepset_candidates: combinations(neigh_y, l)) {
          std::vector<size_t> var_ids(sepset_candidates);
          var_ids.push_back(x);
          var_ids.push_back(y);
          std::sort(var_ids.begin(), var_ids.end());
          ContingencyTable ct = buildContingencyTable(var_ids, df);
          if (citest(x, y, sepset_candidates, ct, citest_type)) {
            pdag.remove_edge(x, y);
            pdag.remove_edge(y, x);
            sepset[x][y] = sepset[y][x] = sepset_candidates;
            break;  // stop checking other sets for this pair
          }
        }
      }
    }
    if (all_nodes_degree_less_than(pdag, l + 1) || l == max_cond_vars) {
        break;
    }
  }
  return pdag;
}

void orient_colliders(
  PDAG &pdag,
  const std::vector<std::vector<std::vector<size_t>>> &sepset) {
  const size_t n = pdag.num_vars;
  for (size_t y = 0; y < n; ++y) {
    auto neigh = pdag.undirected_neighbors(y);
    if (neigh.size() < 2) continue;
    for (size_t i = 0; i < neigh.size(); ++i) {
      for (size_t j = i + 1; j < neigh.size(); ++j) {
        size_t x = neigh[i], z = neigh[j];
        if (pdag.is_adjacent(x, z)) continue;
        if (std::find(sepset[x][z].begin(), sepset[x][z].end(), y) ==
            sepset[x][z].end()) {
          if (pdag.has_edge(y, x)) pdag.remove_edge(y, x);
          if (pdag.has_edge(y, z)) pdag.remove_edge(y, z);
        }
      }
    }
  }
}

PDAG PC(
  const DataframeWrapper &df, 
  const CITestType &citest_type,
  size_t max_cond_vars
) {
  const size_t n = df.num_of_vars;
  std::vector<std::vector<std::vector<size_t>>> sepset(
      n, std::vector<std::vector<size_t>>(n));

  PDAG pdag = build_skeleton(df, citest_type, max_cond_vars, sepset);
  orient_colliders(pdag, sepset);
  return pdag;
}
