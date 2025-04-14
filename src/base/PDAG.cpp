#include "base/PDAG.h"

#include <stdexcept>

PDAG::PDAG(size_t num_vars) : num_vars(num_vars) {
  size_t blocks = (num_vars + 63) / 64;
  adj_mat.resize(num_vars, std::vector<uint64_t>(blocks, 0ULL));
}

PDAG &PDAG::operator=(const PDAG &a) {
  if (this != &a) {
    this->num_vars = a.num_vars;
    this->adj_mat = a.adj_mat;
  }
  return *this;
}

PDAG::PDAG(const PDAG &old) : num_vars(old.num_vars), adj_mat(old.adj_mat) {}

size_t PDAG::get_num_vars() const { return num_vars; }

// if x -> y
bool PDAG::has_edge(size_t x, size_t y) const {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");
  size_t block = y / 64;
  size_t shift = y % 64;
  return (adj_mat[x][block] & (1ULL << shift)) != 0;
}

// if x -> y and not y -> x
bool PDAG::has_directed_edge(size_t x, size_t y) const {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");
  return has_edge(x, y) && !has_edge(y, x);
}

// if x <-> y
bool PDAG::has_undirected_edge(size_t x, size_t y) const {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");
  return has_edge(x, y) && has_edge(y, x);
}

// succ(x) = {y | x -> y}
std::vector<size_t> PDAG::successors(size_t x) const {
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
std::vector<size_t> PDAG::predecessors(size_t x) const {
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
std::vector<size_t> PDAG::neighbors(size_t x) const {
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
std::vector<size_t> PDAG::undirected_neighbors(size_t x) const {
  if (x >= num_vars) throw std::out_of_range("Index out of range");
  std::vector<size_t> undir;
  auto succ = successors(x);
  for (size_t node : succ) {
    if (has_edge(node, x)) undir.push_back(node);
  }
  return undir;
}

bool PDAG::has_directed_path(size_t x, size_t y) const {
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

bool PDAG::has_path(size_t x, size_t y) const {
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

bool PDAG::has_connection(size_t x, size_t y) const {
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
void PDAG::add_edge(size_t x, size_t y) {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");
  if (has_edge(x, y)) throw std::invalid_argument("Edge already exists");
  size_t block = y / 64;
  size_t shift = y % 64;
  adj_mat[x][block] |= (1ULL << shift);
}

// remove edge x -> y
void PDAG::remove_edge(size_t x, size_t y) {
  if (x >= num_vars || y >= num_vars)
    throw std::out_of_range("Index out of range");
  if (!has_edge(x, y)) throw std::invalid_argument("Edge does not exist");
  size_t block = y / 64;
  size_t shift = y % 64;
  adj_mat[x][block] &= ~(1ULL << shift);
}

void PDAG::complete_graph() {
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