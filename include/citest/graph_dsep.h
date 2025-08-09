#pragma once

#include <algorithm>
#include <queue>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "graph/pdag.h"

inline bool is_d_separated(const PDAG& g,
                           std::size_t x,
                           std::size_t y,
                           const std::vector<std::size_t>& Z) {

  // (0) Check input validity
  const std::size_t n = g.num_vars;
  if (x >= n || y >= n) throw std::out_of_range("x or y is out of range");
  for (auto z : Z)
    if (z >= n) throw std::out_of_range("Z contains out-of-range index");
  if (std::find(Z.begin(), Z.end(), x) != Z.end() ||
      std::find(Z.begin(), Z.end(), y) != Z.end()) {
    throw std::invalid_argument("x or y is in Z");
  }
  
  std::vector<char> in_Z(n, 0);
  for (auto z : Z) in_Z[z] = 1;

  // (1) Compute ancestor set An(S) where S = {x, y} ∪ Z
  const auto& parents_of = g.parents;
  std::vector<char> in_anc(n, 0);
  std::queue<std::size_t> q;
  auto push_if_new = [&](std::size_t v) {
    if (!in_anc[v]) {
      in_anc[v] = 1;
      q.push(v);
    }
  };
  push_if_new(x);
  push_if_new(y);
  for (auto z : Z) push_if_new(z);
  while (!q.empty()) {
    auto v = q.front();
    q.pop();
    for (auto p : parents_of[v]) {
      if (!in_anc[p]) {
        in_anc[p] = 1;
        q.push(p);
      }
    }
  }

  // (2) Moralization 
  std::vector<std::unordered_set<std::size_t>> und_adj(n);

  // (2.a) Make directed edges undirected
  for (std::size_t v = 0; v < n; ++v) {
    if (!in_anc[v]) continue;
    for (auto u : parents_of[v]) {
      if (!in_anc[u]) continue;
      und_adj[u].insert(v);
      und_adj[v].insert(u);
    }
  }

  // (2.b) Connect parents of the same child (moral edges)
  for (std::size_t v = 0; v < n; ++v) {
    if (!in_anc[v]) continue;
    std::vector<std::size_t> P;
    P.reserve(parents_of[v].size());
    for (auto p : parents_of[v])
      if (in_anc[p]) P.push_back(p);

    if (P.size() >= 2) {
      for (std::size_t i = 0; i + 1 < P.size(); ++i) {
        for (std::size_t j = i + 1; j < P.size(); ++j) {
          auto a = P[i], b = P[j];
          und_adj[a].insert(b);
          und_adj[b].insert(a);
        }
      }
    }
  }

  // (3) Check reachability from X to Y in the undirected graph (without going through Z)
  std::vector<char> vis(n, 0);
  std::queue<std::size_t> bfs;
  vis[x] = 1;
  bfs.push(x);
  while (!bfs.empty()) {
    auto u = bfs.front();
    bfs.pop();
    for (auto v : und_adj[u]) {
      if (in_Z[v] || vis[v]) continue;
      if (v == y) {
        return false;
      }
      vis[v] = 1;
      bfs.push(v);
    }
  }
  return true; // If we cannot reach y from x, they are d-separated.
}