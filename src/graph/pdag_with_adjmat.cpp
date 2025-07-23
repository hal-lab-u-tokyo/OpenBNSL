#include "graph/pdag_with_adjmat.h"

#include <stdexcept>
#include <unordered_map>

#include "graph/pdag.h"

PDAGwithAdjMat::PDAGwithAdjMat(size_t num_vars) : num_vars(num_vars) {
  size_t blocks = (num_vars + 63) / 64;
  adj_mat.resize(num_vars, std::vector<uint64_t>(blocks, 0ULL));
}

PDAGwithAdjMat &PDAGwithAdjMat::operator=(const PDAGwithAdjMat &a) {
  if (this != &a) {
    this->num_vars = a.num_vars;
    this->adj_mat = a.adj_mat;
  }
  return *this;
}

PDAGwithAdjMat::PDAGwithAdjMat(const PDAGwithAdjMat &old)
    : num_vars(old.num_vars), adj_mat(old.adj_mat) {}

size_t PDAGwithAdjMat::get_num_vars() const { return num_vars; }

// if x -> y
bool PDAGwithAdjMat::has_edge(size_t x, size_t y) const {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");
  size_t block = y / 64;
  size_t shift = y % 64;
  return (adj_mat[x][block] & (1ULL << shift)) != 0;
}

// if x -> y and not y -> x
bool PDAGwithAdjMat::has_directed_edge(size_t x, size_t y) const {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");
  return has_edge(x, y) && !has_edge(y, x);
}

// if x -> y and y -> x
bool PDAGwithAdjMat::has_undirected_edge(size_t x, size_t y) const {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");
  return has_edge(x, y) && has_edge(y, x);
}

// if x -> y or y -> x
bool PDAGwithAdjMat::is_adjacent(size_t x, size_t y) const {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");
  return has_edge(x, y) || has_edge(y, x);
}

// succ(x) = {y | x -> y}
std::vector<size_t> PDAGwithAdjMat::successors(size_t x) const {
  if (x >= num_vars) throw std::out_of_range("Index out of range");
  std::vector<size_t> succ;
  size_t blocks = (num_vars + 63) / 64;
  for (size_t j = 0; j < blocks; ++j) {
    uint64_t bits = adj_mat[x][j];
    while (bits) {
      size_t shift = __builtin_ctzll(bits);
      size_t neighbor = j * 64 + shift;
      if (neighbor >= num_vars) break;
      succ.push_back(neighbor);
      bits &= bits - 1;
    }
  }
  return succ;
}

// pred(x) = {y | y -> x} (slower than successors)
std::vector<size_t> PDAGwithAdjMat::predecessors(size_t x) const {
  if (x >= num_vars) throw std::out_of_range("Index out of range");
  std::vector<size_t> pred;
  size_t block = x / 64;
  size_t shift = x % 64;
  uint64_t mask = (1ULL << shift);
  for (size_t i = 0; i < num_vars; ++i) {
    if (adj_mat[i][block] & mask) pred.push_back(i);
  }
  return pred;
}

// neigh(x) = {y | x -> y or y -> x}
std::vector<size_t> PDAGwithAdjMat::neighbors(size_t x) const {
  if (x >= num_vars) throw std::out_of_range("Index out of range");
  std::vector<size_t> neigh;
  std::vector<bool> visited(num_vars, false);

  auto succ = successors(x);
  for (size_t node : succ) {
    neigh.push_back(node);
    visited[node] = true;
  }
  auto pred = predecessors(x);
  for (size_t node : pred) {
    if (!visited[node]) neigh.push_back(node);
  }
  return neigh;
}

// undirected_neighbors(x) = {y | x <-> y}
// O( |N(x)| )
std::vector<size_t> PDAGwithAdjMat::undirected_neighbors(size_t x) const {
  if (x >= num_vars) throw std::out_of_range("Index out of range");
  std::vector<size_t> undir;
  auto succ = successors(x);
  for (size_t node : succ) {
    if (has_edge(node, x)) undir.push_back(node);
  }
  return undir;
}

bool PDAGwithAdjMat::has_directed_path(size_t x, size_t y) const {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");

  std::vector<bool> visited(num_vars, false);
  std::vector<size_t> stack;
  stack.push_back(x);

  while (!stack.empty()) {
    size_t node = stack.back();
    stack.pop_back();

    if (node == y) return true;
    if (visited[node]) continue;
    visited[node] = true;

    auto succ = successors(node);
    for (auto s : succ) {
      if (!visited[s] && has_directed_edge(node, s)) stack.push_back(s);
    }
  }
  return false;
}

bool PDAGwithAdjMat::has_path(size_t x, size_t y) const {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");

  std::vector<bool> visited(num_vars, false);
  std::vector<size_t> stack;
  stack.push_back(x);

  while (!stack.empty()) {
    size_t node = stack.back();
    stack.pop_back();

    if (node == y) return true;
    if (visited[node]) continue;
    visited[node] = true;

    auto succ = successors(node);
    for (auto s : succ) {
      if (!visited[s]) stack.push_back(s);
    }
  }
  return false;
}

bool PDAGwithAdjMat::has_connection(size_t x, size_t y) const {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");

  std::vector<bool> visited(num_vars, false);
  std::vector<size_t> stack;
  stack.push_back(x);

  while (!stack.empty()) {
    size_t node = stack.back();
    stack.pop_back();

    if (node == y) return true;
    if (visited[node]) continue;
    visited[node] = true;

    auto neigh = neighbors(node);
    for (auto n : neigh) {
      if (!visited[n]) stack.push_back(n);
    }
  }
  return false;
}

// add edge x -> y
void PDAGwithAdjMat::add_edge(size_t x, size_t y) {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");
  if (has_edge(x, y)) throw std::invalid_argument("Edge already exists");
  size_t block = y / 64;
  size_t shift = y % 64;
  adj_mat[x][block] |= (1ULL << shift);
}

// remove edge x -> y
void PDAGwithAdjMat::remove_edge(size_t x, size_t y) {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");
  if (!has_edge(x, y)) throw std::invalid_argument("Edge does not exist");
  size_t block = y / 64;
  size_t shift = y % 64;
  adj_mat[x][block] &= ~(1ULL << shift);
}

void PDAGwithAdjMat::complete_graph() {
  size_t blocks = (num_vars + 63) / 64;
  for (size_t i = 0; i < num_vars; ++i) {
    for (size_t j = 0; j < blocks; ++j) {
      if (j == blocks - 1) {
        if (num_vars % 64 == 0)
          adj_mat[i][j] = ~0ULL;
        else
          adj_mat[i][j] = (1ULL << (num_vars % 64)) - 1;
      } else {
        adj_mat[i][j] = ~0ULL;
      }
    }
  }
  for (size_t i = 0; i < num_vars; ++i) {
    remove_edge(i, i);
  }
}

PDAG PDAGwithAdjMat::to_pdag() const {
  PDAG pdag(num_vars);
  for (size_t i = 0; i < num_vars; ++i) {
    for (size_t j = 0; j < num_vars; ++j) {
      if (has_edge(i, j)) {
        pdag.add_edge(i, j);
      }
    }
  }
  return pdag;
}

// TODO
// void apply_meeks_rules(PDAGwithAdjMat &g) {
//   const size_t n = g.num_vars;
//   bool changed = true;

//   while (changed) {
//     changed = false;

//     // ---------------------------- Rule 1 ----------------------------
//     for (size_t y = 0; y < n; ++y) {
//       auto parents = directed_parents(g, y);
//       auto undirected_nbs = g.undirected_neighbors(g, y);
//       for (size_t x : parents) {
//         for (size_t z : undirected_nbs) {
//           if (is_adjacent(g, x, z)) continue;
//           // Avoid cycles: ensure there is no directed path z ⇒ y
//           if (has_directed_path(g, z, y)) continue;
//           // Orient y‑z as y→z
//           if (g.has_edge(z, y)) {
//             g.remove_edge(z, y);
//             changed = true;
//           }
//         }
//       }
//     }

//     // ---------------------------- Rule 2 ----------------------------
//     for (size_t z = 0; z < n; ++z) {
//       auto parents = directed_parents(g, z);
//       auto children = directed_children(g, z);
//       for (size_t x : parents) {
//         for (size_t y : children) {
//           if (g.has_undirected_edge(g, x, y)) {
//             // Orient x‑y as x→y
//             g.remove_edge(y, x);
//             changed = true;
//           }
//         }
//       }
//     }

//     // ---------------------------- Rule 3 ----------------------------
//     for (size_t x = 0; x < n; ++x) {
//       auto undirected_nbs = g.undirected_neighbors(g, x);
//       if (undirected_nbs.size() < 3) continue;
//       for (size_t w : undirected_nbs) {
//         size_t cnt = 0;
//         for (size_t nbr : undirected_nbs) {
//           if (nbr == w) continue;
//           if (g.has_directed_edge(g, nbr, w)) ++cnt;
//           if (cnt >= 2) break;
//         }
//         if (cnt >= 2) {
//           // Orient x‑w as x→w
//           if (g.has_edge(w, x)) {
//             g.remove_edge(w, x);
//             changed = true;
//           }
//         }
//       }
//     }
//   }
// }
