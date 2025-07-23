#pragma once

#include <cstdint>
#include <vector>

#include "graph/ipdag_convertible.h"
#include "graph/pdag.h"

struct PDAGwithAdjMat : IPDAGConvertible {
  std::size_t num_vars;
  std::vector<std::vector<uint64_t>> adj_mat;

  PDAGwithAdjMat(std::size_t num_vars);
  PDAGwithAdjMat(const PDAGwithAdjMat &old);
  PDAGwithAdjMat &operator=(const PDAGwithAdjMat &a);
  ~PDAGwithAdjMat() = default;

  std::size_t get_num_vars() const;
  std::vector<std::size_t> successors(std::size_t x) const;    // {y | x -> y}
  std::vector<std::size_t> predecessors(std::size_t x) const;  // {y | y -> x}
  std::vector<std::size_t> neighbors(
      std::size_t x) const;  // {y | x -> y or y -> x}
  std::vector<std::size_t> undirected_neighbors(
      std::size_t x) const;                               // {y | x <-> y}
  bool has_edge(std::size_t from, std::size_t to) const;  // x -> y
  bool has_directed_edge(std::size_t from,
                         std::size_t to) const;  // x -> y and not y -> x
  bool has_undirected_edge(std::size_t from,
                           std::size_t to) const;  // x -> y and y -> x
  bool is_adjacent(std::size_t from, std::size_t to) const;  // x -> y or y -> x
  bool has_directed_path(std::size_t from, std::size_t to) const;
  bool has_path(std::size_t from, std::size_t to) const;
  bool has_connection(std::size_t from, std::size_t to) const;
  void add_edge(std::size_t from, std::size_t to);
  void remove_edge(std::size_t from, std::size_t to);
  void complete_graph();

  PDAG to_pdag() const;
};