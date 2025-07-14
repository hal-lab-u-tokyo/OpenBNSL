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
  std::vector<size_t> successors(size_t x) const;    // {y | x -> y}
  std::vector<size_t> predecessors(size_t x) const;  // {y | y -> x}
  std::vector<size_t> neighbors(size_t x) const;     // {y | x -> y or y -> x}
  std::vector<size_t> undirected_neighbors(size_t x) const;  // {y | x <-> y}
  bool has_edge(size_t from, size_t to) const; // x -> y
  bool has_directed_edge(size_t from, size_t to) const; // x -> y and not y -> x
  bool has_undirected_edge(size_t from, size_t to) const; // x -> y and y -> x
  bool is_adjacent(size_t from, size_t to) const;  // x -> y or y -> x
  bool has_directed_path(size_t from, size_t to) const;
  bool has_path(size_t from, size_t to) const;
  bool has_connection(size_t from, size_t to) const;
  void add_edge(size_t from, size_t to);
  void remove_edge(size_t from, size_t to);
  void complete_graph();
};