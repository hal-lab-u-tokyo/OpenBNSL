#pragma once
#include <cstddef>
#include <set>
#include <unordered_set>
#include <vector>

#include "graph/ipdag_convertible.h"
#include "graph/pdag.h"

template <bool Deterministic>
using ParentSet = std::conditional_t<Deterministic,
                                     std::set<std::size_t>,
                                     std::unordered_set<std::size_t>>;

/**
 * @ingroup graph
 * @struct PDAGwithAdjList
 * @brief PDAG implementation backed by an adjacency list.
 * @details
 * This structure uses an adjacency list to represent the PDAG, allowing for
 * efficient storage and access patterns.
 */
template <bool Deterministic>
struct PDAGwithAdjList : IPDAGConvertible {
  using ParentSetType = ParentSet<Deterministic>;

  std::size_t num_vars;
  std::vector<ParentSet<Deterministic>> parents;

  explicit PDAGwithAdjList(std::size_t num_vars);
  bool has_edge(std::size_t from, std::size_t to) const;
  void add_edge(std::size_t from, std::size_t to);
  void remove_edge(std::size_t from, std::size_t to);
  bool has_path(std::size_t src, std::size_t dst) const;
  void set_parents(std::size_t v, const ParentSetType& new_parents);
  PDAG to_pdag() const;
};