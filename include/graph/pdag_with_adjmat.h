#pragma once

#include <cstdint>
#include <limits>
#include <unordered_set>
#include <vector>

#include "graph/ipdag_convertible.h"
#include "graph/pdag.h"
#include "utils/logging.h"

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
  std::size_t num_vars; // < 1000
  std::vector<std::vector<uint64_t>> potential_parent_bits;
  
  PDAGwithAdjMat(std::size_t n);
  void set_as_complete();

  PDAG to_pdag() const override;
  void orient_colliders(const Sepset& sepset);
  void apply_meeks_rules();

  bool has_directed_edge(std::size_t u, std::size_t v) const {
    return _has_arc(u, v) && !_has_arc(v, u);
  }
  bool has_undirected_edge(std::size_t u, std::size_t v) const {
    return _has_arc(u, v) && _has_arc(v, u);
  }
  bool is_adjacent(std::size_t u, std::size_t v) const {
    return _has_arc(u, v) || _has_arc(v, u);
  }
  void orient_edge(std::size_t u, std::size_t v) {
    if (!is_adjacent(u, v)) throw std::invalid_argument("No edge exists");
    if (has_directed_edge(v, u)) 
      // throw std::invalid_argument("Edge already oriented the other way");
      INFO("[PDAG] edge already oriented the other way");
    _set_arc(u, v);
    _clr_arc(v, u);
  }
  void remove_edge(std::size_t u, std::size_t v) {
    _clr_arc(u, v);
    _clr_arc(v, u);
  }

  /*
   * potential_parents(v) = {u | u->v or u-v}
   * parents(v) = {u | u->v}
   * undirected_neighbors(v) = {u | u-v}
   */
  std::vector<std::size_t> potential_parents(std::size_t v) const;
  std::vector<std::size_t> parents(std::size_t v) const;
  std::vector<std::size_t> undirected_neighbors(std::size_t v) const;
  std::vector<std::size_t> potential_parents_in(
    std::size_t v, const std::vector<uint64_t>& mask) const;
  std::vector<std::size_t> parents_in(
    std::size_t v, const std::vector<uint64_t>& mask) const;
  std::vector<std::size_t> undirected_neighbors_in(
    std::size_t v, const std::vector<uint64_t>& mask) const;

  inline bool _has_arc(std::size_t u, std::size_t v) const {
    const std::size_t b = u / 64, s = u % 64;
    return (potential_parent_bits[v][b] & (1ULL << s)) != 0ULL;
  }
  inline void _set_arc(std::size_t u, std::size_t v) {
    const std::size_t b = u / 64, s = u % 64;
    potential_parent_bits[v][b] |= (1ULL << s);
  }
  inline void _clr_arc(std::size_t u, std::size_t v) {
    const std::size_t b = u / 64, s = u % 64;
    potential_parent_bits[v][b] &= ~(1ULL << s);
  }
  inline std::vector<size_t>
  _extract_indices_from_bits(const std::vector<uint64_t>& bits) const {
    std::vector<size_t> out;
    const std::size_t blocks = (num_vars + 63) / 64;
    for (std::size_t b = 0; b < blocks; ++b) {
      uint64_t x = bits[b];
      while (x) {
  #if defined(__GNUC__) || defined(__clang__)
        std::size_t s = __builtin_ctzll(x);
  #else
        unsigned long s; _BitScanForward64(&s, x);
  #endif
        const std::size_t id = b * 64 + s;
        if (id < num_vars) out.push_back(id);
        x &= x - 1;
      }
    }
    return out;
  }
};