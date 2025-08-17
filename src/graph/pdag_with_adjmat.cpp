#include "graph/pdag_with_adjmat.h"

#include <algorithm>
#include <stdexcept>

/* Lifecycle */
PDAGwithAdjMat::PDAGwithAdjMat(size_t num_vars) : num_vars(num_vars) {
  size_t blocks = (num_vars + 63) / 64;
  adj_mat.resize(num_vars, std::vector<uint64_t>(blocks, 0ULL));
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
bool PDAGwithAdjMat::_has_arc(std::size_t u, std::size_t v) const {
  _bounds(u, num_vars);
  _bounds(v, num_vars);
  std::size_t b = v / 64, s = v % 64;
  return (adj_mat[u][b] & (1ULL << s)) != 0ULL;
}
void PDAGwithAdjMat::_remove_arc(std::size_t u, std::size_t v) {
  _bounds(u, num_vars);
  _bounds(v, num_vars);
  std::size_t b = v / 64, s = v % 64;
  adj_mat[u][b] &= ~(1ULL << s);
}

/* Read-only primitives */
bool PDAGwithAdjMat::has_directed_edge(std::size_t u, std::size_t v) const {
  return _has_arc(u, v) && !_has_arc(v, u);
}
bool PDAGwithAdjMat::has_undirected_edge(std::size_t u, std::size_t v) const {
  return _has_arc(u, v) && _has_arc(v, u);
}
bool PDAGwithAdjMat::is_adjacent(std::size_t u, std::size_t v) const {
  return _has_arc(u, v) || _has_arc(v, u);
}

/* Modification primitives */
void PDAGwithAdjMat::remove_undirected_edge(std::size_t u, std::size_t v) {
  if (!has_undirected_edge(u, v)) throw std::invalid_argument("edge not found");
  _remove_arc(u, v);
  _remove_arc(v, u);
}
void PDAGwithAdjMat::orient_edge(std::size_t u, std::size_t v) {
  if (has_directed_edge(u, v)) return;  // already oriented
  _remove_arc(v, u);
}

/* Neighbor operations */
// {u | v -> u or v <-> u}
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
  return res;  // {u | v -> u or v <-> u}
}
std::vector<std::size_t> PDAGwithAdjMat::children(std::size_t v) const {
  _bounds(v, num_vars);
  std::vector<std::size_t> res;
  for (auto u : successors(v))  // v -> u or v <-> u
    if (!_has_arc(u, v)) res.push_back(u);
  return res;  // {u | v -> u}
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors(
    std::size_t v) const {
  std::vector<std::size_t> res;
  for (auto u : successors(v))  // v -> u or v <-> u
    if (_has_arc(u, v)) res.push_back(u);
  return res;  // {u | v <-> u}
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors_without(
    std::size_t v,
    std::size_t excl) const {
  std::vector<std::size_t> res;
  for (auto u : successors(v))  // v -> u or v <-> u
    if (u != excl && _has_arc(u, v)) res.push_back(u);
  return res;  // {u | v <-> u and u != excl}
}

/* Reachability operations */

/* Graph-wide operations */
PDAG PDAGwithAdjMat::to_pdag() const {
  PDAG g(num_vars);
  for (std::size_t u = 0; u < num_vars; ++u)
    for (std::size_t v = 0; v < num_vars; ++v)
      if (_has_arc(u, v)) g.add_edge(u, v);
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
    _remove_arc(i, i);  // no self-loop
  }
}

/* Meek's rules and helpers */
void PDAGwithAdjMat::apply_meeks_rules(bool apply_r4) {
  bool changed = true;
  while (changed) {
    changed = false;

    /* ---------- Rule 1 ---------- */
    // If (1) A->B, (2) B-C, and (3) A and C are non-adjacent, then orient B->C.
    // Reason: otherwise B<-C would induce an unexpected unshielded collider
    // A->B<-C.
    for (std::size_t a = 0; a < num_vars; ++a) {
      for (std::size_t b : children(a)) {
        for (std::size_t c : undirected_neighbors_without(b, a)) {
          if (!is_adjacent(a, c)) {
            orient_edge(b, c);
            changed = true;
          }
        }
      }
    }

    /* ---------- Rule 2 ---------- */
    // If (1) A->B, (2) B->C, and (3) A-C, then orient A->C.
    // Reason: otherwise A<-C would create a directed cycle A->B->C->A.
    for (std::size_t a = 0; a < num_vars; ++a) {
      for (std::size_t b : children(a)) {
        for (std::size_t c : children(b)) {
          if (has_undirected_edge(a, c)) {
            orient_edge(a, c);
            changed = true;
          }
        }
      }
    }

    /* ---------- Rule 3 ---------- */
    // If (1) A-{B,C,D}, (2) {B,C} -> D and (3) B and C are non-adjacent, then
    // orient A->D. Reason: if we instead chose D->A, then to avoid cycles we
    // would be forced to orient B->A and C->A. But since B and C are
    // non-adjacent, this would introduce a new unshielded collider B->A<-C.
    for (std::size_t a = 0; a < num_vars; ++a) {
      auto bcd_candidates = undirected_neighbors(a);
      if (bcd_candidates.size() < 3) continue;  // require at least 3 neighbors

      // For each unordered pair (b,c) of parents with b and c non-adjacent
      for (std::size_t j = 0; j < bcd_candidates.size(); ++j) {
        for (std::size_t k = j + 1; k < bcd_candidates.size(); ++k) {
          std::size_t b = bcd_candidates[j];
          std::size_t c = bcd_candidates[k];
          if (is_adjacent(b, c)) continue;
          for (std::size_t d : bcd_candidates) {
            if (d == b || d == c) continue;
            if (has_directed_edge(b, d) && has_directed_edge(c, d)) {
              orient_edge(a, d);
              changed = true;
            }
          }
        }
      }
    }

    /* ---------- Rule 4 ---------- */
    // Out of scope
  }
}