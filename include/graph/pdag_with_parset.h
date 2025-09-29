#pragma once

#include <cstddef>
#include <limits>
#include <unordered_set>
#include <vector>

#include "graph/ipdag_convertible.h"
#include "graph/pdag.h"

constexpr size_t UNASSIGNED = std::numeric_limits<size_t>::max();

using Sepset = std::vector<std::vector<std::unordered_set<size_t>>>;

/**
 * @ingroup graph
 * @struct PDAGwithParSet
 */
struct PDAGwithParSet : IPDAGConvertible {
  std::size_t num_local_vars;
  std::size_t num_global_vars;
  std::vector<size_t> gid2lid;  // global id -> local idx or UNASSIGNED
  std::vector<size_t> lid2gid;  // local idx -> global id

  std::vector<std::unordered_set<size_t>> in_parents;

  explicit PDAGwithParSet(std::size_t n);
  static PDAGwithParSet induced_subgraph(const PDAGwithParSet& src,
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

  /*
   * predecessors(v) = {u | u->v or u<->v}
   * parents(v) = {u | u->v}
   * undirected_neighbors(v) = {u | u-v}
   * undirected_neighbors_without(v, excl) = {u | u-v and u!=excl}
   */
  std::unordered_set<std::size_t> predecessors(std::size_t vG) const;
  std::unordered_set<std::size_t> parents(std::size_t vG) const;
  std::unordered_set<std::size_t> undirected_neighbors(std::size_t vG) const;
  std::unordered_set<std::size_t> undirected_neighbors_without(
      std::size_t vG,
      std::size_t exclG) const;

  std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>> decompose(
      const PDAGwithParSet& g_all) const;

  void orient_colliders(const Sepset& sepset);
  void apply_meeks_rules();

  PDAG to_pdag() const override;

 private:
  static inline void _bounds(std::size_t i, std::size_t n) {
    if (i >= n) throw std::out_of_range("PDAGwithParSet index");
  }
  inline std::size_t _to_local(std::size_t g) const {
    return (g < gid2lid.size()) ? gid2lid[g] : UNASSIGNED;
  }

  bool has_directed_edgeL(std::size_t uL, std::size_t vL) const {
    return in_parents[vL].count(uL) && !in_parents[uL].count(vL);
  }
  bool has_undirected_edgeL(std::size_t uL, std::size_t vL) const {
    return in_parents[vL].count(uL) && in_parents[uL].count(vL);
  }
  bool is_adjacentL(std::size_t uL, std::size_t vL) const {
    return in_parents[vL].count(uL) || in_parents[uL].count(vL);
  }
  void remove_undirected_edgeL(std::size_t uL, std::size_t vL) {
    if (!has_undirected_edgeL(uL, vL)) return;
    in_parents[vL].erase(uL);
    in_parents[uL].erase(vL);
  }
  void orient_edgeL(std::size_t uL, std::size_t vL) {
    // if (!is_adjacentL(uL, vL)) return;
    // in_parents[vL].insert(uL);
    // in_parents[uL].erase(vL);
    if (has_directed_edgeL(uL, vL)) return;
    in_parents[uL].erase(vL);
  }

  std::unordered_set<std::size_t> predecessorsL(std::size_t vL) const;
  std::unordered_set<std::size_t> parentsL(std::size_t vL) const;
  std::unordered_set<std::size_t> undirected_neighborsL(std::size_t vL) const;
  std::unordered_set<std::size_t> undirected_neighbors_withoutL(
      std::size_t vL,
      std::size_t exclL) const;
};
