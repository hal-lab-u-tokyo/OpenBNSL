#include "graph/pdag_with_adjmat.h"

#include <algorithm>
#include <queue>
#include <stdexcept>
#include "utils/logging.h"

PDAGwithAdjMat::PDAGwithAdjMat(std::size_t num_vars) : num_vars(num_vars) {
  const std::size_t blocks = (num_vars + 63) / 64;
  potential_parent_bits.assign(num_vars, std::vector<uint64_t>(blocks, 0ULL));
}

void PDAGwithAdjMat::set_as_complete() {
  const std::size_t blocks = (num_vars + 63) / 64;
  for (std::size_t i = 0; i < num_vars; ++i) {
    for (std::size_t b = 0; b < blocks; ++b) {
      const bool tail = (b == blocks - 1) && (num_vars % 64);
      const uint64_t mask =
          tail ? ((1ULL << (num_vars % 64)) - 1) : ~0ULL;
      potential_parent_bits[i][b] = mask;
    }
    _clr_arc(i, i);  // no self-loop
  }
}

PDAG PDAGwithAdjMat::to_pdag() const {
  PDAG p(num_vars);
  for (std::size_t v = 0; v < num_vars; ++v) {
    const auto p_pars = potential_parents(v);
    for (auto u : p_pars) {
      p.add_edge(u, v);
      if (_has_arc(v, u)) p.add_edge(v, u);
    }
  }
  return p;
}

std::vector<std::size_t> PDAGwithAdjMat::potential_parents(std::size_t v) const {
  return _extract_indices_from_bits(potential_parent_bits[v]);
}

std::vector<std::size_t> PDAGwithAdjMat::parents(std::size_t v) const {
  std::vector<std::size_t> res, src = potential_parents(v);
  for (auto u : src) {
    if (!_has_arc(v, u)) res.push_back(u);
  }
  return res;
}

std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors(std::size_t v) const {
  std::vector<std::size_t> res, src = potential_parents(v);
  for (auto u : src) {
    if (_has_arc(v, u)) res.push_back(u);
  }
  return res;
}

std::vector<std::size_t> PDAGwithAdjMat::potential_parents_in(
    std::size_t v, const std::vector<uint64_t>& mask) const {
  const std::size_t blocks = (num_vars + 63) / 64;
  std::vector<uint64_t> bits(blocks);
  for (std::size_t b = 0; b < blocks; ++b) {
    bits[b] = potential_parent_bits[v][b] & mask[b];
  }
  return _extract_indices_from_bits(bits);
}

std::vector<std::size_t> PDAGwithAdjMat::parents_in(
    std::size_t v, const std::vector<uint64_t>& mask) const {
  std::vector<std::size_t> res;
  const auto src = potential_parents_in(v, mask);
  for (auto u : src) {
    if (!_has_arc(v, u)) res.push_back(u);
  }
  return res;
}

std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors_in(
    std::size_t v, const std::vector<uint64_t>& mask) const {
  std::vector<std::size_t> res;
  const auto src = potential_parents_in(v, mask);
  for (auto u : src) {
    if (_has_arc(v, u)) res.push_back(u);
  }
  return res;
}


/* ---------- Rule V ----------
 * If (1) y-{x,z}, (2) x and z are non-adjacent, and
 * (3) y is not in the sepset of x and z, then orient x->y<-z.
 */
void PDAGwithAdjMat::orient_colliders(const Sepset& sepset) {
  std::unordered_set<uint64_t> orient_pairs;
  for (std::size_t y = 0; y < num_vars; ++y) {
    const auto p_pars = potential_parents(y); // {x | x->y or x-y}
    for (std::size_t i = 0; i + 1 < p_pars.size(); ++i) {
      for (std::size_t j = i + 1; j < p_pars.size(); ++j) {
        const auto x = p_pars[i], z = p_pars[j];
        if (is_adjacent(x, z)) continue;  // shielded
        if (sepset[x][z].count(y) == 0) {
          uint64_t pair1 = (static_cast<uint64_t>(x) << 32) | (static_cast<uint64_t>(y));
          uint64_t pair2 = (static_cast<uint64_t>(z) << 32) | (static_cast<uint64_t>(y));
          orient_pairs.insert(pair1);
          orient_pairs.insert(pair2);
        }
      }
    }
  }
  for (auto pair : orient_pairs) {
    auto p = static_cast<std::size_t>(pair >> 32);
    auto c = static_cast<std::size_t>(pair & 0xFFFFFFFF);
    orient_edge(p, c);
  }
  INFO("[PDAG] oriented " << orient_pairs.size() << " colliders");
}

void PDAGwithAdjMat::apply_meeks_rules() {
restart:
  /* ---------- Rule 1 ----------
    * If (1) x->y, (2) y-z, and (3) x and z are non-adjacent, then orient y->z.
    */
  for (std::size_t y = 0; y < num_vars; ++y) {
    const auto x_candidates = parents(y);
    const auto z_candidates = undirected_neighbors(y);
    for (auto x : x_candidates) {
      for (auto z : z_candidates) {
        if (!is_adjacent(x, z)) {
          orient_edge(y, z); // y->z
          goto restart;
        }
      }
    }
  }

  /* ---------- Rule 2 ----------
    * If (1) x->y, (2) y->z, and (3) x-z, then orient x->z.
    */
  for (std::size_t z = 0; z < num_vars; ++z) {
    const auto x_candidates = undirected_neighbors(z);
    const auto y_candidates = parents(z);
    for (auto x : x_candidates) {
      for (auto y : y_candidates) {
        if (has_directed_edge(x, y)) {
          orient_edge(x, z);
          goto restart;
        }
      }
    }
  }

  /* ---------- Rule 3 ----------
    * If (1) x-{y,z,w}, (2) {y,z} -> w, and (3) y and z are non-adjacent,
    * then orient x->w.
    */
  for (std::size_t w = 0; w < num_vars; ++w) {
    const auto x_candidates = undirected_neighbors(w);
    const auto yz_candidates = parents(w);
    for (std::size_t i = 0; i + 1 < yz_candidates.size(); ++i) {
      for (std::size_t j = i + 1; j < yz_candidates.size(); ++j) {
        const auto y = yz_candidates[i], z = yz_candidates[j];
        if (is_adjacent(y, z)) continue;
        for (auto x : x_candidates) {
          if (has_undirected_edge(x, y) && has_undirected_edge(x, z)) {
            orient_edge(x, w);
            goto restart;
          }
        }
      }
    }
  }
}
