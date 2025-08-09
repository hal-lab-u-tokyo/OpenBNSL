#include "graph/pdag_with_adjmat.h"

#include <algorithm>
#include <stdexcept>

/* Lifecycle */
PDAGwithAdjMat::PDAGwithAdjMat(size_t num_vars) : num_vars(num_vars) {
  size_t blocks = (num_vars + 63) / 64;
  adj_mat.resize(num_vars, std::vector<uint64_t>(blocks, 0ULL));
  complete_graph();
}
PDAGwithAdjMat::PDAGwithAdjMat(const PDAGwithAdjMat &old)
    : num_vars(old.num_vars), adj_mat(old.adj_mat) {}
PDAGwithAdjMat &PDAGwithAdjMat::operator=(const PDAGwithAdjMat &a) {
  if (this != &a) {
    this->num_vars = a.num_vars;
    this->adj_mat = a.adj_mat;
  }
  return *this;
}

/* Private helpers */
inline static void _bounds(std::size_t idx, std::size_t n) {
  if (idx >= n) throw std::out_of_range("PDAG index out of range");
}

/* Read-only primitives */
bool PDAGwithAdjMat::has_directed_edge(std::size_t u, std::size_t v) const {
  _bounds(u, num_vars);
  _bounds(v, num_vars);
  std::size_t b = v / 64, s = v % 64;
  return (adj_mat[u][b] & (1ULL << s)) != 0ULL;
}
bool PDAGwithAdjMat::has_undirected_edge(std::size_t u, std::size_t v) const {
  return has_directed_edge(u, v) && has_directed_edge(v, u);
}
bool PDAGwithAdjMat::is_adjacent(std::size_t u, std::size_t v) const {
  return has_directed_edge(u, v) || has_directed_edge(v, u);
}

/* Modification primitives */
void PDAGwithAdjMat::add_directed_edge(std::size_t u, std::size_t v) {
  _bounds(u, num_vars);
  _bounds(v, num_vars);
  if (has_directed_edge(u, v))
    throw std::invalid_argument("arc already exists");
  adj_mat[u][v / 64] |= (1ULL << (v % 64));
}
void PDAGwithAdjMat::remove_directed_edge(std::size_t u, std::size_t v) {
  _bounds(u, num_vars);
  _bounds(v, num_vars);
  if (!has_directed_edge(u, v)) throw std::invalid_argument("arc not found");
  adj_mat[u][v / 64] &= ~(1ULL << (v % 64));
}
void PDAGwithAdjMat::add_undirected_edge(std::size_t u, std::size_t v) {
  if (has_undirected_edge(u, v))
    throw std::invalid_argument("edge already exists");
  add_directed_edge(u, v);
  add_directed_edge(v, u);
}
void PDAGwithAdjMat::remove_undirected_edge(std::size_t u, std::size_t v) {
  if (!has_undirected_edge(u, v)) throw std::invalid_argument("edge not found");
  remove_directed_edge(u, v);
  remove_directed_edge(v, u);
}
void PDAGwithAdjMat::orient_edge(std::size_t from, std::size_t to) {
  if (!has_undirected_edge(from, to)) return;
  remove_directed_edge(to, from);
}

/* Neighbor operations */
std::vector<std::size_t> PDAGwithAdjMat::successors(std::size_t v) const {
  _bounds(v, num_vars);
  std::vector<std::size_t> res;
  std::size_t blocks = (num_vars + 63) / 64;
  for (std::size_t j = 0; j < blocks; ++j) {
    uint64_t bits = adj_mat[v][j];
    while (bits) {
      std::size_t s = __builtin_ctzll(bits);
      std::size_t idx = j * 64 + s;
      if (idx < num_vars) res.push_back(idx);
      bits &= bits - 1;
    }
  }
  return res;
}
std::vector<std::size_t> PDAGwithAdjMat::predecessors(std::size_t v) const {
  _bounds(v, num_vars);
  std::vector<std::size_t> res;
  std::size_t b = v / 64, m = 1ULL << (v % 64);
  for (std::size_t u = 0; u < num_vars; ++u)
    if (adj_mat[u][b] & m) res.push_back(u);
  return res;
}
std::vector<std::size_t> PDAGwithAdjMat::neighbors(std::size_t v) const {
  auto res = successors(v);
  auto pred = predecessors(v);
  res.insert(res.end(), pred.begin(), pred.end());
  std::sort(res.begin(), res.end());
  res.erase(std::unique(res.begin(), res.end()), res.end());
  return res;
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors(
    std::size_t v) const {
  std::vector<std::size_t> res;
  for (auto u : successors(v))
    if (has_undirected_edge(u, v)) res.push_back(u);
  return res;
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors_without(
    std::size_t v,
    std::size_t excl) const {
  std::vector<std::size_t> res;
  for (auto u : successors(v))
    if (u != excl && has_undirected_edge(u, v)) res.push_back(u);
  return res;
}

/* Reachability operations */
bool PDAGwithAdjMat::has_directed_path(std::size_t s, std::size_t t) const {
  _bounds(s, num_vars);
  _bounds(t, num_vars);
  std::vector<char> vis(num_vars);
  std::vector<std::size_t> st{s};
  while (!st.empty()) {
    auto u = st.back();
    st.pop_back();
    if (u == t) return true;
    if (vis[u]) continue;
    vis[u] = 1;
    for (auto w : successors(u))
      if (!vis[w] && !has_directed_edge(w, u))
        st.push_back(w);  // only forward arcs
  }
  return false;
}
bool PDAGwithAdjMat::has_path(std::size_t s, std::size_t t) const {
  _bounds(s, num_vars);
  _bounds(t, num_vars);
  std::vector<char> vis(num_vars);
  std::vector<std::size_t> st{s};
  while (!st.empty()) {
    auto u = st.back();
    st.pop_back();
    if (u == t) return true;
    if (vis[u]) continue;
    vis[u] = 1;
    for (auto w : successors(u))
      if (!vis[w]) st.push_back(w);
  }
  return false;
}
bool PDAGwithAdjMat::has_connection(std::size_t s, std::size_t t) const {
  _bounds(s, num_vars);
  _bounds(t, num_vars);
  std::vector<char> vis(num_vars);
  std::vector<std::size_t> st{s};
  while (!st.empty()) {
    auto u = st.back();
    st.pop_back();
    if (u == t) return true;
    if (vis[u]) continue;
    vis[u] = 1;
    for (auto w : neighbors(u))
      if (!vis[w]) st.push_back(w);
  }
  return false;
}

/* Graph-wide operations */
PDAG PDAGwithAdjMat::to_pdag() const {
  PDAG g(num_vars);
  for (std::size_t u = 0; u < num_vars; ++u)
    for (std::size_t v = 0; v < num_vars; ++v)
      if (has_directed_edge(u, v)) g.add_edge(u, v);
  return g;
}
void PDAGwithAdjMat::complete_graph() {
  std::size_t blocks = (num_vars + 63) / 64;
  for (std::size_t i = 0; i < num_vars; ++i) {
    for (std::size_t j = 0; j < blocks; ++j) {
      uint64_t mask = (j == blocks - 1 && num_vars % 64)
                          ? ((1ULL << (num_vars % 64)) - 1)
                          : ~0ULL;
      adj_mat[i][j] = mask;
    }
    remove_directed_edge(i, i);  // no self-loop
  }
}

/* Meek's rules and helpers */
std::vector<std::size_t> PDAGwithAdjMat::directed_parents(std::size_t v) const {
  std::vector<std::size_t> res;
  for (auto u : predecessors(v))
    if (!has_directed_edge(v, u)) res.push_back(u);  // u → v かつ v ↛ u
  return res;
}
std::vector<std::size_t> PDAGwithAdjMat::directed_children(
    std::size_t v) const {
  std::vector<std::size_t> res;
  for (auto u : successors(v))
    if (!has_directed_edge(u, v)) res.push_back(u);  // v → u かつ u ↛ v
  return res;
}
std::vector<std::size_t> PDAGwithAdjMat::all_neighbors(std::size_t v) const {
  auto res = neighbors(v);
  return res;
}
bool PDAGwithAdjMat::creates_unshielded_collider(std::size_t y,
                                                 std::size_t z) const {
  for (auto x : directed_parents(y))
    if (!is_adjacent(x, z)) return true;
  return false;
}
void PDAGwithAdjMat::apply_meeks_rules(bool apply_r4) {
  bool changed = true;
  while (changed) {
    changed = false;

    /* ---------- Rule 1 ---------- */
    for (std::size_t y = 0; y < num_vars; ++y)
      for (auto x : directed_parents(y))
        for (auto z : undirected_neighbors(y))
          if (!is_adjacent(x, z) && !creates_unshielded_collider(y, z) &&
              !has_directed_path(z, y)) {
            orient_edge(y, z);
            changed = true;
          }

    /* ---------- Rule 2 ---------- */
    for (std::size_t z = 0; z < num_vars; ++z)
      for (auto x : directed_parents(z))
        for (auto y : directed_children(z))
          if (has_undirected_edge(x, y)) {
            orient_edge(x, y);
            changed = true;
          }

    /* ---------- Rule 3 ---------- */
    for (std::size_t x = 0; x < num_vars; ++x) {
      auto und = undirected_neighbors(x);
      if (und.size() < 3) continue;

      for (auto y : und)
        for (auto z : und)
          if (y != z)
            for (auto w : und)
              if (w != y && w != z)
                if (has_directed_edge(y, w) && !has_directed_edge(w, y) &&
                    has_directed_edge(z, w) && !has_directed_edge(w, z)) {
                  orient_edge(x, w);
                  changed = true;
                  goto NEXT_X;
                }
    NEXT_X:;
    }

    /* ---------- Rule 4 (optional) ---------- */
    if (apply_r4) {
      for (std::size_t c = 0; c < num_vars; ++c)
        for (auto b : directed_children(c))
          for (auto d : directed_parents(c)) {
            if (b == d || is_adjacent(b, d)) continue;
            std::vector<std::size_t> cand = undirected_neighbors(b);
            auto tmp = undirected_neighbors(d);
            cand.insert(cand.end(), tmp.begin(), tmp.end());
            tmp = all_neighbors(c);
            cand.insert(cand.end(), tmp.begin(), tmp.end());
            std::sort(cand.begin(), cand.end());
            cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
            for (auto a : cand)
              if (has_undirected_edge(a, b)) {
                orient_edge(a, b);
                changed = true;
              }
          }
    }
  }
}