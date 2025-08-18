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
  std::size_t b = u / 64, s = u % 64;
  return (adj_mat[v][b] & (1ULL << s)) != 0ULL;
}
void PDAGwithAdjMat::_remove_arc(std::size_t u, std::size_t v) {
  _bounds(u, num_vars);
  _bounds(v, num_vars);
  std::size_t b = u / 64, s = u % 64;
  adj_mat[v][b] &= ~(1ULL << s);
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
void PDAGwithAdjMat::remove_directed_edge(std::size_t u, std::size_t v) {
  if (!has_directed_edge(u, v)) throw std::invalid_argument("edge not found");
  _remove_arc(u, v);
}
void PDAGwithAdjMat::orient_edge(std::size_t u, std::size_t v) {
  if (has_directed_edge(u, v)) return;  // already oriented
  _remove_arc(v, u);
}

/* Neighbor operations */
std::vector<std::size_t> PDAGwithAdjMat::predecessors(std::size_t v) const {
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
  return res;  // {u | v <- u or v <-> u}
}

std::vector<std::size_t> PDAGwithAdjMat::parents(std::size_t v) const {
  _bounds(v, num_vars);
  std::vector<std::size_t> res;
  for (auto u : predecessors(v))  // v <- u or v <-> u
    if (!_has_arc(v, u)) res.push_back(u);
  return res;  // {u | v <- u}
}

std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors(
    std::size_t v) const {
  std::vector<std::size_t> res;
  for (auto u : predecessors(v))  // v <- u or v <-> u
    if (_has_arc(v, u)) res.push_back(u);
  return res;  // {u | v <-> u}
}

std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors_without(
    std::size_t v,
    std::size_t excl) const {
  std::vector<std::size_t> res;
  for (auto u : predecessors(v))  // v <- u or v <-> u
    if (u != excl && _has_arc(v, u)) res.push_back(u);
  return res;  // {u | v <-> u and u != excl}
}

/* Reachability operations */

/* Graph-wide operations */
PDAG PDAGwithAdjMat::to_pdag() const {
  PDAG g(num_vars);
  for (std::size_t u = 0; u < num_vars; ++u)
    for (std::size_t v = 0; v < num_vars; ++v)
      if (_has_arc(v, u)) g.add_edge(v, u);
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

void PDAGwithAdjMat::orient_colliders(const Sepset &sepset) {
  for (size_t y = 0; y < num_vars; ++y) {
    bool changed = true;
    while (changed) {
      changed = false;
      auto neigh = undirected_neighbors(y);
      for (size_t i = 0; i + 1 < neigh.size(); ++i) {
        for (size_t j = i + 1; j < neigh.size(); ++j) {
          size_t x = neigh[i], z = neigh[j];
          if (is_adjacent(x, z)) continue;
          if (sepset[x][z].count(y) == 0) {
            orient_edge(x, y);
            orient_edge(z, y);
            changed = true;
          }
        }
      }
    }
  }
}

void PDAGwithAdjMat::apply_meeks_rules() {
  bool changed = true;
  while (changed) {
    changed = false;

    /* ---------- Rule 1 ---------- */
    // If (1) A->B, (2) B-C, and (3) A and C are non-adjacent, then orient B->C.
    // Reason: orient B<-C would induce an unexpected unshielded collider
    // A->B<-C.
    for (std::size_t b = 0; b < num_vars; ++b) {
      auto a_candidates = parents(b);               // { a | a -> b }
      auto c_candidates = undirected_neighbors(b);  // { c | b <-> c }
      for (std::size_t a : a_candidates) {
        for (std::size_t c : c_candidates) {
          if (!is_adjacent(a, c)) {
            orient_edge(b, c);
            changed = true;
          }
        }
      }
    }

    /* ---------- Rule 2 ---------- */
    // If (1) A->B, (2) B->C, and (3) A-C, then orient A->C.
    // Reason: oriente A<-C would create a directed cycle A->B->C->A.
    for (std::size_t c = 0; c < num_vars; ++c) {
      auto a_candidates = undirected_neighbors(c);  // { a | a <-> c }
      auto b_candidates = parents(c);               // { b | b -> c }
      for (std::size_t a : a_candidates) {
        for (std::size_t b : b_candidates) {
          if (has_directed_edge(a, b)) {
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
    for (std::size_t d = 0; d < num_vars; ++d) {
      auto a_candidates = undirected_neighbors(d);
      auto bc_candidates = parents(d);
      for (std::size_t i = 0; i < bc_candidates.size(); ++i) {
        for (std::size_t j = i + 1; j < bc_candidates.size(); ++j) {
          std::size_t b = bc_candidates[i];
          std::size_t c = bc_candidates[j];
          if (is_adjacent(b, c)) continue;
          for (std::size_t a : a_candidates) {
            if (has_undirected_edge(a, b) && has_undirected_edge(a, c)) {
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