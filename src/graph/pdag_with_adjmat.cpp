#include "graph/pdag_with_adjmat.h"

#include <algorithm>
#include <queue>
#include <stdexcept>

inline void PDAGwithAdjMat::_bounds(std::size_t i, std::size_t n) {
  if (i >= n) throw std::out_of_range("PDAGwithAdjMat index");
}

PDAGwithAdjMat::PDAGwithAdjMat(std::size_t n)
    : num_local_vars(n), num_global_vars(n) {
  gid2lid.assign(num_global_vars, UNASSIGNED);
  lid2gid.assign(num_local_vars, UNASSIGNED);

  const std::size_t blocks = (num_local_vars + 63) / 64;
  in_arc_bits.assign(num_local_vars, std::vector<uint64_t>(blocks, 0ULL));

  for (std::size_t i = 0; i < num_local_vars; ++i) {
    gid2lid[i] = i;
    lid2gid[i] = i;
  }
}

PDAGwithAdjMat PDAGwithAdjMat::induced_subgraph(const PDAGwithAdjMat& src,
                                                const std::vector<size_t>& S) {
  PDAGwithAdjMat res(S.size());
  res.num_global_vars = src.num_global_vars;
  res.gid2lid.assign(res.num_global_vars, UNASSIGNED);
  res.lid2gid.assign(res.num_local_vars, UNASSIGNED);
  for (std::size_t i = 0; i < res.num_local_vars; ++i) {
    const size_t g = S[i];
    res.gid2lid[g] = i;
    res.lid2gid[i] = g;
  }

  const std::size_t blocks = (res.num_local_vars + 63) / 64;
  res.in_arc_bits.assign(res.num_local_vars,
                         std::vector<uint64_t>(blocks, 0ULL));
  for (std::size_t i = 0; i < res.num_local_vars; ++i) {
    const size_t g_i = res.lid2gid[i];
    for (auto g_nb : src.undirected_neighbors(g_i)) {
      const size_t j = res._to_local(g_nb);
      if (j == UNASSIGNED) continue;
      res._set_arcL(i, j);
      res._set_arcL(j, i);
    }
  }
  return res;
}

bool PDAGwithAdjMat::_has_arcL(std::size_t uL, std::size_t vL) const {
  _bounds(uL, num_local_vars);
  _bounds(vL, num_local_vars);
  const std::size_t b = uL / 64, s = uL % 64;
  return (in_arc_bits[vL][b] & (1ULL << s)) != 0ULL;
}
void PDAGwithAdjMat::_set_arcL(std::size_t uL, std::size_t vL) {
  _bounds(uL, num_local_vars);
  _bounds(vL, num_local_vars);
  const std::size_t b = uL / 64, s = uL % 64;
  in_arc_bits[vL][b] |= (1ULL << s);
}
void PDAGwithAdjMat::_clr_arcL(std::size_t uL, std::size_t vL) {
  _bounds(uL, num_local_vars);
  _bounds(vL, num_local_vars);
  const std::size_t b = uL / 64, s = uL % 64;
  in_arc_bits[vL][b] &= ~(1ULL << s);
}

void PDAGwithAdjMat::set_as_complete() {
  const std::size_t blocks = (num_local_vars + 63) / 64;
  for (std::size_t i = 0; i < num_local_vars; ++i) {
    for (std::size_t b = 0; b < blocks; ++b) {
      const bool tail = (b == blocks - 1) && (num_local_vars % 64);
      const uint64_t mask =
          tail ? ((1ULL << (num_local_vars % 64)) - 1) : ~0ULL;
      in_arc_bits[i][b] = mask;
    }
    _clr_arcL(i, i);  // no self-loop
  }
}

PDAG PDAGwithAdjMat::to_pdag() const {
  PDAG p(num_global_vars);
  for (std::size_t vL = 0; vL < num_local_vars; ++vL) {
    for (auto uL : predecessorsL(vL)) {
      p.add_edge(lid2gid[uL], lid2gid[vL]);
      if (_has_arcL(vL, uL)) p.add_edge(lid2gid[vL], lid2gid[uL]);
    }
  }
  return p;
}

std::vector<std::size_t> PDAGwithAdjMat::predecessorsL(std::size_t vL) const {
  _bounds(vL, num_local_vars);
  std::vector<std::size_t> out;
  const std::size_t blocks = (num_local_vars + 63) / 64;
  for (std::size_t b = 0; b < blocks; ++b) {
    uint64_t bits = in_arc_bits[vL][b];
    while (bits) {
#if defined(__GNUC__) || defined(__clang__)
      std::size_t s = __builtin_ctzll(bits);
#else
      unsigned long idx;
      _BitScanForward64(&idx, bits);  // <intrin.h>
      std::size_t s = static_cast<std::size_t>(idx);
#endif
      std::size_t uL = b * 64 + s;
      if (uL < num_local_vars) out.push_back(uL);
      bits &= bits - 1;  // clear LSB
    }
  }
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::parentsL(std::size_t vL) const {
  std::vector<std::size_t> out;
  for (auto uL : predecessorsL(vL)) {
    if (!_has_arcL(vL, uL)) out.push_back(uL);
  }
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighborsL(
    std::size_t vL) const {
  std::vector<std::size_t> out;
  for (auto uL : predecessorsL(vL)) {
    if (_has_arcL(vL, uL)) out.push_back(uL);
  }
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors_withoutL(
    std::size_t vL,
    std::size_t excl) const {
  std::vector<std::size_t> out;
  for (auto uL : predecessorsL(vL)) {
    if (uL != excl && _has_arcL(vL, uL)) out.push_back(uL);
  }
  return out;
}

bool PDAGwithAdjMat::has_directed_edge(std::size_t uG, std::size_t vG) const {
  const auto uL = _to_local(uG), vL = _to_local(vG);
  if (uL == UNASSIGNED || vL == UNASSIGNED) return false;
  return has_directed_edgeL(uL, vL);
}
bool PDAGwithAdjMat::has_undirected_edge(std::size_t uG, std::size_t vG) const {
  const auto uL = _to_local(uG), vL = _to_local(vG);
  if (uL == UNASSIGNED || vL == UNASSIGNED) return false;
  return has_undirected_edgeL(uL, vL);
}
bool PDAGwithAdjMat::is_adjacent(std::size_t uG, std::size_t vG) const {
  const auto uL = _to_local(uG), vL = _to_local(vG);
  if (uL == UNASSIGNED || vL == UNASSIGNED) return false;
  return is_adjacentL(uL, vL);
}

void PDAGwithAdjMat::remove_undirected_edge(std::size_t uG, std::size_t vG) {
  const auto uL = _to_local(uG), vL = _to_local(vG);
  if (uL == UNASSIGNED || vL == UNASSIGNED) return;  // not in this subgraph
  remove_undirected_edgeL(uL, vL);
}
void PDAGwithAdjMat::orient_edge(std::size_t uG, std::size_t vG) {
  const auto uL = _to_local(uG), vL = _to_local(vG);
  if (uL == UNASSIGNED || vL == UNASSIGNED) return;
  orient_edgeL(uL, vL);
}

std::vector<std::size_t> PDAGwithAdjMat::predecessors(std::size_t vG) const {
  const auto vL = _to_local(vG);
  if (vL == UNASSIGNED) throw std::runtime_error("v not in this graph");
  std::vector<size_t> out;
  for (auto uL : predecessorsL(vL)) out.push_back(lid2gid[uL]);
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::parents(std::size_t vG) const {
  const auto vL = _to_local(vG);
  if (vL == UNASSIGNED) throw std::runtime_error("v not in this graph");
  std::vector<size_t> out;
  for (auto uL : parentsL(vL)) out.push_back(lid2gid[uL]);
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors(
    std::size_t vG) const {
  const auto vL = _to_local(vG);
  if (vL == UNASSIGNED) throw std::runtime_error("v not in this graph");
  std::vector<size_t> out;
  for (auto uL : undirected_neighborsL(vL)) out.push_back(lid2gid[uL]);
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors_without(
    std::size_t vG,
    std::size_t exclG) const {
  const auto vL = _to_local(vG);
  if (vL == UNASSIGNED) throw std::runtime_error("v not in this graph");
  const auto exclL = _to_local(exclG);
  std::vector<size_t> out;
  for (auto uL : undirected_neighbors_withoutL(vL, exclL))
    out.push_back(lid2gid[uL]);
  return out;
}

std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>>
PDAGwithAdjMat::decompose(const PDAGwithAdjMat& g_all) const {
  std::vector<size_t> desc_nodes_local;
  std::vector<size_t> desc_nodes;
  std::vector<size_t> remaining_nodes_local;

  for (std::size_t xL = 0; xL < num_local_vars; ++xL) {
    bool has_child = false;
    for (std::size_t yL = 0; yL < num_local_vars; ++yL) {
      if (has_directed_edgeL(xL, yL)) {
        has_child = true;
        break;
      }
    }
    if (has_child) {
      remaining_nodes_local.push_back(xL);
    } else {
      desc_nodes_local.push_back(xL);
      desc_nodes.push_back(lid2gid[xL]);
    }
  }

  std::vector<std::vector<size_t>> asc_nodes_list;
  std::vector<char> seen(num_local_vars, 0);

  for (auto sL : remaining_nodes_local) {
    if (seen[sL]) continue;

    std::vector<size_t> asc_nodes;
    std::queue<size_t> q;
    q.push(sL);
    seen[sL] = 1;

    while (!q.empty()) {
      const auto uL = q.front();
      q.pop();
      const auto uG = lid2gid[uL];
      asc_nodes.push_back(uG);

      for (auto vL : remaining_nodes_local) {
        if (seen[vL]) continue;
        const auto vG = lid2gid[vL];
        if (g_all.is_adjacent(uG, vG)) {
          seen[vL] = 1;
          q.push(vL);
        }
      }
    }
    asc_nodes_list.push_back(std::move(asc_nodes));
  }

  return {desc_nodes, asc_nodes_list};
}

/* ---------- Rule V ---------- */
/*
 * If (1) y-{x,z}, (2) x and z are non-adjacent, and
 * (3) y is not in the sepset of x and z, then orient x->y<-z.
 */
void PDAGwithAdjMat::orient_colliders(const Sepset& sepset) {
  for (std::size_t yL = 0; yL < num_local_vars; ++yL) {
    const auto pres = predecessorsL(yL);  // {u | u->y or u<->y}
    for (std::size_t i = 0; i + 1 < pres.size(); ++i) {
      for (std::size_t j = i + 1; j < pres.size(); ++j) {
        const auto xL = pres[i], zL = pres[j];
        if (is_adjacentL(xL, zL)) continue;  // shielded
        const size_t xG = lid2gid[xL], yG = lid2gid[yL], zG = lid2gid[zL];
        if (sepset[xG][zG].count(yG) == 0) {
          orient_edgeL(xL, yL);
          orient_edgeL(zL, yL);
        }
      }
    }
  }
}

void PDAGwithAdjMat::apply_meeks_rules() {
  bool changed = true;
  while (changed) {
    changed = false;

    /* ---------- Rule 1 ----------
     * If (1) x->y, (2) y-z, and (3) x and z are non-adjacent, then orient y->z.
     */
    for (std::size_t yL = 0; yL < num_local_vars; ++yL) {
      const auto xL_candidates = parentsL(yL);
      const auto zL_candidates = undirected_neighborsL(yL);
      for (auto xL : xL_candidates) {
        for (auto zL : zL_candidates) {
          if (!is_adjacentL(xL, zL)) {
            orient_edgeL(yL, zL);
            changed = true;
          }
        }
      }
    }

    /* ---------- Rule 2 ----------
     * If (1) x->y, (2) y->z, and (3) x-z, then orient x->z.
     */
    for (std::size_t zL = 0; zL < num_local_vars; ++zL) {
      const auto xL_candidates = undirected_neighborsL(zL);
      const auto yL_candidates = parentsL(zL);
      for (auto xL : xL_candidates) {
        for (auto yL : yL_candidates) {
          if (has_directed_edgeL(xL, yL)) {
            orient_edgeL(xL, zL);
            changed = true;
          }
        }
      }
    }

    /* ---------- Rule 3 ----------
     * If (1) x-{y,z,w}, (2) {y,z} -> w, and (3) y and z are non-adjacent,
     * then orient x->w.
     */
    for (std::size_t wL = 0; wL < num_local_vars; ++wL) {
      const auto xL_candidates = undirected_neighborsL(wL);
      const auto yzL_candidates = parentsL(wL);
      for (std::size_t i = 0; i + 1 < yzL_candidates.size(); ++i) {
        for (std::size_t j = i + 1; j < yzL_candidates.size(); ++j) {
          const auto yL = yzL_candidates[i], zL = yzL_candidates[j];
          if (is_adjacentL(yL, zL)) continue;
          for (auto xL : xL_candidates) {
            if (has_undirected_edgeL(xL, yL) && has_undirected_edgeL(xL, zL)) {
              orient_edgeL(xL, wL);
              changed = true;
            }
          }
        }
      }
    }
  }
}
