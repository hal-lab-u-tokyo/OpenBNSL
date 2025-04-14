#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

struct PDAG {
  size_t num_vars;
  std::vector<std::vector<uint64_t>> adj_mat;

  PDAG(size_t num_vars);
  PDAG(const PDAG &old);
  PDAG &operator=(const PDAG &a);
  ~PDAG() = default;

  size_t get_num_vars() const;
  bool has_edge(size_t from, size_t to) const;
  void add_edge(size_t from, size_t to);
  void remove_edge(size_t from, size_t to);
  void complete_graph();
  std::vector<size_t> neighbors(size_t x) const;
};