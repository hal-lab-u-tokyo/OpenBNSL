#pragma once

#include <algorithm>
#include <queue>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "graph/pdag.h"

// PDAG について：
// - G.has_edge(u, v) は u -> v の有向辺があるかを返します。
// - 無向辺は (u->v) と (v->u) の両方が存在するものとして表現されます。
// - d-separation 判定は DAG
// 前提のアルゴリズム（祖先グラフのモラライズ）を用います。
//   CPDAG
//   を与えた場合は、両向き辺がある箇所は無向として扱われ、標準のモラライズ基準に一致します。

inline bool is_d_separated(const PDAG& G,
                           std::size_t X,
                           std::size_t Y,
                           const std::vector<std::size_t>& Z) {
  const std::size_t n = G.num_vars;
  if (X >= n || Y >= n) throw std::out_of_range("X or Y is out of range");
  for (auto z : Z)
    if (z >= n) throw std::out_of_range("Z contains out-of-range index");

  if (X == Y) {
    // 同一点なら「分離されていない」と解釈して false を返す（活性経路が
    // trivially ある）
    return false;
  }

  // Z をフラグ化
  std::vector<char> in_Z(n, 0);
  for (auto z : Z) in_Z[z] = 1;

  // 1) 祖先集合 An(S) の計算（S = {X, Y} ∪ Z）
  //    parents[v] = {u | u -> v} を上向きに辿る
  const auto& parents_of = G.parents;

  std::vector<char> in_anc(n, 0);
  std::queue<std::size_t> q;
  auto push_if_new = [&](std::size_t v) {
    if (!in_anc[v]) {
      in_anc[v] = 1;
      q.push(v);
    }
  };

  push_if_new(X);
  push_if_new(Y);
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

  // 2) モラライズ：祖先グラフ上で無向隣接表を作る
  std::vector<std::unordered_set<std::size_t>> und_adj(n);

  // 2a) 親→子の有向辺を無向化（u - v）
  for (std::size_t v = 0; v < n; ++v) {
    if (!in_anc[v]) continue;
    for (auto u : parents_of[v]) {
      if (!in_anc[u]) continue;
      und_adj[u].insert(v);
      und_adj[v].insert(u);
    }
  }

  // 2b) 同一子の親同士を接続（moral edges）
  for (std::size_t v = 0; v < n; ++v) {
    if (!in_anc[v]) continue;
    // 祖先グラフ内の親だけ抽出
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

  // 3) Z を削除（探索で訪問禁止）
  if (in_Z[X] || in_Z[Y]) {
    // 片方でも条件付けされていれば、端点でブロックされるので分離
    return true;
  }

  // 4) 無向グラフで X→Y の到達可否（Z を通らない）
  std::vector<char> vis(n, 0);
  std::queue<std::size_t> bfs;
  vis[X] = 1;
  bfs.push(X);

  while (!bfs.empty()) {
    auto u = bfs.front();
    bfs.pop();
    for (auto v : und_adj[u]) {
      if (in_Z[v] || vis[v]) continue;
      if (v == Y) {
        // 到達できる → d-connected（非分離）
        return false;
      }
      vis[v] = 1;
      bfs.push(v);
    }
  }

  // 到達できなければ d-separated
  return true;
}