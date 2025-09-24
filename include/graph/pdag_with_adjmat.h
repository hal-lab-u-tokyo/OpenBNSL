#pragma once

#include <cstdint>
#include <limits>
#include <unordered_set>
#include <vector>

#include "graph/ipdag_convertible.h"
#include "graph/pdag.h"

constexpr size_t UNASSIGNED = std::numeric_limits<size_t>::max();

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
  std::size_t num_local_vars;
  std::size_t num_global_vars;
  std::vector<size_t> gid2lid;  // global id -> local idx or UNASSIGNED
  std::vector<size_t> lid2gid;  // local idx -> global id

  std::vector<std::vector<uint64_t>> in_arc_bits;

  explicit PDAGwithAdjMat(std::size_t n);
  static PDAGwithAdjMat induced_subgraph(const PDAGwithAdjMat& src,
                                         const std::vector<size_t>& S);

  void set_as_complete();
  bool contains(size_t gid) const {
    return gid < gid2lid.size() && gid2lid[gid] != UNASSIGNED;
  }

  bool has_directed_edge(std::size_t uG, std::size_t vG) const;
  bool has_undirected_edge(std::size_t uG, std::size_t vG) const;
  bool is_adjacent(std::size_t uG, std::size_t vG) const;

  void remove_undirected_edge(std::size_t uG, std::size_t vG);
  void orient_edge(std::size_t uG, std::size_t vG);

  std::vector<std::size_t> predecessors(
      std::size_t vG) const;  // {uG | u->v or u<->v}
  std::vector<std::size_t> parents(std::size_t vG) const;  // {uG | u->v}
  std::vector<std::size_t> undirected_neighbors(
      std::size_t vG) const;  // {uG | u-v}
  std::vector<std::size_t> undirected_neighbors_without(
      std::size_t vG,
      std::size_t exclG) const;

  std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>> decompose(
      const PDAGwithAdjMat& g_all) const;

  void orient_colliders(const Sepset& sepset);
  void apply_meeks_rules();

  PDAG to_pdag() const override;

 private:
  static inline void _bounds(std::size_t i, std::size_t n);

  inline std::size_t _to_local(std::size_t g) const {
    return (g < gid2lid.size()) ? gid2lid[g] : UNASSIGNED;
  }

  bool _has_arcL(std::size_t uL, std::size_t vL) const;
  void _set_arcL(std::size_t uL, std::size_t vL);
  void _clr_arcL(std::size_t uL, std::size_t vL);

  bool has_directed_edgeL(std::size_t uL, std::size_t vL) const {
    return _has_arcL(uL, vL) && !_has_arcL(vL, uL);
  }
  bool has_undirected_edgeL(std::size_t uL, std::size_t vL) const {
    return _has_arcL(uL, vL) && _has_arcL(vL, uL);
  }
  bool is_adjacentL(std::size_t uL, std::size_t vL) const {
    return _has_arcL(uL, vL) || _has_arcL(vL, uL);
  }
  void remove_undirected_edgeL(std::size_t uL, std::size_t vL) {
    if (!has_undirected_edgeL(uL, vL)) return;
    _clr_arcL(uL, vL);
    _clr_arcL(vL, uL);
  }
  void orient_edgeL(std::size_t uL, std::size_t vL) {  // u -> v にする
    if (has_directed_edgeL(uL, vL)) return;
    _clr_arcL(vL, uL);
  }

  std::vector<std::size_t> predecessorsL(std::size_t vL) const;
  std::vector<std::size_t> parentsL(std::size_t vL) const;
  std::vector<std::size_t> undirected_neighborsL(std::size_t vL) const;
  std::vector<std::size_t> undirected_neighbors_withoutL(
      std::size_t vL,
      std::size_t excl) const;
};