#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

class PDAG {
 private:
  size_t n;
  std::vector<std::vector<uint64_t>> adj_mat;

 public:
  PDAG(size_t n);
  PDAG(const PDAG &old);
  PDAG &operator=(const PDAG &a);
  ~PDAG() = default;

  bool has_edge(size_t from, size_t to) const;
  void add_edge(size_t from, size_t to);
  void remove_edge(size_t from, size_t to);
  void complete_graph();
  std::vector<size_t> neighbors(size_t x) const;
};