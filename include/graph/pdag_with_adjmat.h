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
  std::size_t num_vars;
  std::vector<size_t> g2l_map;  // global id -> local idx or UNASSIGNED
  std::vector<size_t> var_ids;  // local idx -> global id
  std::vector<std::vector<uint64_t>> adj_mat;

  PDAGwithAdjMat(std::size_t n);
  static PDAGwithAdjMat induced_subgraph(const PDAGwithAdjMat& G,
                                         const std::vector<size_t>& S);
  void set_as_complete();

  bool has_arc(std::size_t u, std::size_t v) const;
  void set_arc(std::size_t u, std::size_t v);
  void clr_arc(std::size_t u, std::size_t v);

  bool has_directed_edge(std::size_t u, std::size_t v) const {
    return has_arc(u, v) && !has_arc(v, u);
  }
  bool has_undirected_edge(std::size_t u, std::size_t v) const {
    return has_arc(u, v) && has_arc(v, u);
  }
  bool is_adjacent(std::size_t u, std::size_t v) const {
    return has_arc(u, v) || has_arc(v, u);
  }
  void remove_undirected_edge(std::size_t u, std::size_t v) {
    if (!has_undirected_edge(u, v)) return;  // TODO: error?
    clr_arc(u, v);
    clr_arc(v, u);
  }
  void orient_edge(std::size_t u, std::size_t v) {
    if (has_directed_edge(u, v)) return;  // already oriented
    clr_arc(v, u);
  }

  std::vector<std::size_t> predecessors(std::size_t v) const;
  std::vector<std::size_t> parents(std::size_t v) const;
  std::vector<std::size_t> undirected_neighbors(std::size_t v) const;
  std::vector<std::size_t> undirected_neighbors_without(std::size_t v,
                                                        std::size_t excl) const;

  void orient_colliders(const Sepset& sepset);
  void apply_meeks_rules();

  std::vector<size_t> childless_nodes() const;
  PDAG to_pdag() const override;
};
