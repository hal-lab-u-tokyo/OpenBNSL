#pragma once

#include <cstdint>
#include <unordered_set>
#include <vector>

#include "graph/ipdag_convertible.h"
#include "graph/pdag.h"

using Sepset = std::vector<std::vector<std::unordered_set<size_t>>>;

/**
 * @ingroup graph
 * @struct PDAGwithAdjMat
 * @brief PDAG implementation backed by a bit-compressed adjacency matrix.
 * @details
 * This structure uses a bit-compressed adjacency matrix to represent the PDAG,
 * allowing for efficient storage and access patterns.
 */
struct PDAGwithAdjMat : IPDAGConvertible {
  /* Data members */
  std::size_t num_vars;
  std::vector<std::vector<uint64_t>> adj_mat;

  /* Lifecycle */
  /**
   * @brief Construct a new PDAGwithAdjMat with a specified number of variables.
   * @param num_vars The number of variables in the PDAG.
   */
  PDAGwithAdjMat(std::size_t num_vars);
  PDAGwithAdjMat(const PDAGwithAdjMat &old);
  PDAGwithAdjMat &operator=(const PDAGwithAdjMat &a);
  ~PDAGwithAdjMat() = default;

  /* Private helpers */
  bool _has_arc(std::size_t u, std::size_t v) const;
  void _remove_arc(std::size_t u, std::size_t v);

  /* Read-only operations */
  bool has_directed_edge(std::size_t u, std::size_t v) const;  // u -> v
  bool has_undirected_edge(std::size_t u,
                           std::size_t v) const;         // u -> v and u <- v
  bool is_adjacent(std::size_t u, std::size_t v) const;  // u -> v or  u <- v

  /* Neighbor operations */
  std::vector<std::size_t> predecessors(
      std::size_t v) const;  // {u | v <- u or v <-> u}
  std::vector<std::size_t> parents(std::size_t v) const;  // {u | v <- u}
  std::vector<std::size_t> undirected_neighbors(
      std::size_t v) const;  // {u | v <-> u}
  std::vector<std::size_t> undirected_neighbors_without(
      std::size_t v,
      std::size_t excl) const;  // {u | v <-> u and u != excl}

  /* Reachability operations */

  /* Modification operations */
  void remove_undirected_edge(std::size_t u, std::size_t v);  // delete u <-> v
  void remove_directed_edge(std::size_t u, std::size_t v);    // delete u -> v
  void orient_edge(std::size_t u, std::size_t v);  // (u <-> v) => (u -> v)

  /* Graph-wide operations */
  PDAG to_pdag() const;
  void complete_graph();
  void orient_colliders(const Sepset &sepset);
  void apply_meeks_rules();
};