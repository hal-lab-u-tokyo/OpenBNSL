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
  PDAG to_pdag() const override;

  bool _has_arc(std::size_t uL, std::size_t vL) const;
  void _set_arc(std::size_t uL, std::size_t vL);
  void _clr_arc(std::size_t uL, std::size_t vL);

  bool has_directed_edgeL(std::size_t uL, std::size_t vL) const {
    return _has_arc(uL, vL) && !_has_arc(vL, uL);
  }
  bool has_undirected_edgeL(std::size_t uL, std::size_t vL) const {
    return _has_arc(uL, vL) && _has_arc(vL, uL);
  }
  bool is_adjacentL(std::size_t uL, std::size_t vL) const {
    return _has_arc(uL, vL) || _has_arc(vL, uL);
  }
  void remove_undirected_edgeL(std::size_t uL, std::size_t vL) {
    if (!has_undirected_edgeL(uL, vL)) return;  // TODO: error?
    _clr_arc(uL, vL);
    _clr_arc(vL, uL);
  }
  void orient_edgeL(std::size_t uL, std::size_t vL) {
    if (has_directed_edgeL(uL, vL)) return;  // already oriented
    _clr_arc(vL, uL);
  }

  std::vector<std::size_t> predecessorsL(std::size_t vL) const;
  std::vector<std::size_t> parentsL(std::size_t vL) const;
  std::vector<std::size_t> undirected_neighborsL(std::size_t vL) const;
  std::vector<std::size_t> undirected_neighbors_withoutL(std::size_t vL,
                                                        std::size_t excl) const;

  std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>> decompose(
      const PDAGwithAdjMat& g_all) const;

  void orient_colliders(const Sepset& sepset);
  void apply_meeks_rules();
};
