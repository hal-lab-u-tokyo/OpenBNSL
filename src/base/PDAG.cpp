#include "base/PDAG.h"

#include <stdexcept>

PDAG::PDAG(size_t num_vars) : num_vars(num_vars) {
  size_t blocks = (num_vars + 63) / 64;
  adj_mat.resize(num_vars, std::vector<uint64_t>(blocks, 0ULL));
}

PDAG::PDAG(const PDAG &old) : num_vars(old.num_vars), adj_mat(old.adj_mat) {}

size_t PDAG::get_num_vars() const { return num_vars; }

PDAG &PDAG::operator=(const PDAG &a) {
  if (this != &a) {
    this->num_vars = a.num_vars;
    this->adj_mat = a.adj_mat;
  }
  return *this;
}

bool PDAG::has_edge(size_t from, size_t to) const {
  if (from >= num_vars || to >= num_vars)
    throw std::out_of_range("Index out of range");
  size_t block = to / 64;
  size_t shift = to % 64;
  return (adj_mat[from][block] & (1ULL << shift)) != 0;
}

void PDAG::add_edge(size_t from, size_t to) {
  if (from >= num_vars || to >= num_vars)
    throw std::out_of_range("Index out of range");
  if (has_edge(from, to)) throw std::invalid_argument("Edge already exists");
  size_t block = to / 64;
  size_t shift = to % 64;
  adj_mat[from][block] |= (1ULL << shift);
}

void PDAG::remove_edge(size_t from, size_t to) {
  if (from >= num_vars || to >= num_vars)
    throw std::out_of_range("Index out of range");
  if (!has_edge(from, to)) throw std::invalid_argument("Edge does not exist");
  size_t block = to / 64;
  size_t shift = to % 64;
  adj_mat[from][block] &= ~(1ULL << shift);
}

void PDAG::complete_graph() {
  size_t blocks = (num_vars + 63) / 64;
  for (size_t i = 0; i < num_vars; ++i) {
    for (size_t j = 0; j < blocks; ++j) {
      if (j == blocks - 1) {
        adj_mat[i][j] =
            (1ULL << (num_vars % 64)) - 1;  // Set bits for the last block
      } else {
        adj_mat[i][j] = ~0ULL;  // Set all bits to 1
      }
    }
  }
  for (size_t i = 0; i < num_vars; ++i) {
    remove_edge(i, i);  // Remove self-loops
  }
}

std::vector<size_t> PDAG::neighbors(size_t node) const {
  if (node >= num_vars) throw std::out_of_range("Index out of range");
  size_t blocks = (num_vars + 63) / 64;
  std::vector<size_t> neigh;
  for (size_t j = 0; j < blocks; ++j) {
    uint64_t bits = adj_mat[node][j];
    while (bits) {
      size_t shift =
          __builtin_ctzll(bits);  // Get the index of the least significant bit
      size_t neighbor = j * 64 + shift;
      if (neighbor >= num_vars) throw std::out_of_range("Index out of range");
      neigh.push_back(neighbor);
      bits &= bits - 1;  // Clear the least significant bit
    }
  }
  return neigh;
}