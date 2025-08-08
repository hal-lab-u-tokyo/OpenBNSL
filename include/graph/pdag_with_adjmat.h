#pragma once

#include <cstdint>
#include <vector>

#include "graph/ipdag_convertible.h"
#include "graph/pdag.h"

/**
 * @brief PDAG implementation backed by a bit-compressed adjacency matrix.
 *
 *  *arc*  : directed link  (u -> v)
 *  *edge* : undirected link (u -> v and u <- v) == u <-> v
 */
struct PDAGwithAdjMat : IPDAGConvertible {
  /* Data members */
  std::size_t num_vars;
  std::vector<std::vector<uint64_t>> adj_mat;

  /* Lifecycle */
  PDAGwithAdjMat(std::size_t num_vars);
  PDAGwithAdjMat(const PDAGwithAdjMat &old);
  PDAGwithAdjMat &operator=(const PDAGwithAdjMat &a);
  ~PDAGwithAdjMat() = default;

  /* Graph-wide operations */
  PDAG to_pdag() const;
  void complete_graph();

  /* Read-only operations */
  bool has_directed_edge(std::size_t u, std::size_t v) const;  // u -> v
  bool has_undirected_edge(std::size_t u,
                           std::size_t v) const;         // u -> v and u <- v
  bool is_adjacent(std::size_t u, std::size_t v) const;  // u -> v or  u <- v

  /* Neighbor operations */
  std::vector<std::size_t> successors(std::size_t v) const;    // {w | v -> w}
  std::vector<std::size_t> predecessors(std::size_t v) const;  // {w | w -> v}
  std::vector<std::size_t> neighbors(
      std::size_t v) const;  // {w | v -> w or w -> v}
  std::vector<std::size_t> undirected_neighbors(
      std::size_t v) const;  // {w | v <-> w}
  std::vector<std::size_t> undirected_neighbors_without(
      std::size_t v,
      std::size_t excl) const;  // {w | v <-> w and w != excl}

  /* Reachability operations */
  bool has_directed_path(std::size_t from,
                         std::size_t to) const;           // from -> ... -> to
  bool has_path(std::size_t from, std::size_t to) const;  // ignores orientation
  bool has_connection(std::size_t from, std::size_t to) const;  // ≒ has_path

  /* Modification operations */
  void add_directed_edge(std::size_t u, std::size_t v);       // u -> v
  void add_undirected_edge(std::size_t u, std::size_t v);     // u <-> v
  void remove_directed_edge(std::size_t u, std::size_t v);    // delete u -> v
  void remove_undirected_edge(std::size_t u, std::size_t v);  // delete u <-> v
  void orient_edge(std::size_t from,
                   std::size_t to);  // (from <-> to) → (from -> to)

  /* Meek's rules and helpers */
  std::vector<std::size_t> directed_parents(std::size_t v) const;  // w → v only
  std::vector<std::size_t> directed_children(
      std::size_t v) const;  // v → w only
  std::vector<std::size_t> all_neighbors(
      std::size_t v) const;  // union of pred+succ
  bool creates_unshielded_collider(std::size_t y, std::size_t z) const;
  void apply_meeks_rules(bool apply_r4 = false);
};