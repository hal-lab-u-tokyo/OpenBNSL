#include "graph/pdag_with_adjmat.h"

#include <algorithm>
#include <limits>
#include <queue>
#include <stdexcept>
constexpr size_t UNASSIGNED = std::numeric_limits<size_t>::max();

static inline void _bounds(std::size_t i, std::size_t n) {
  if (i >= n) throw std::out_of_range("PDAGwithAdjMat index");
}

PDAGwithAdjMat::PDAGwithAdjMat(std::size_t n) : num_vars(n) {
  g2l_map.assign(n, UNASSIGNED);
  var_ids.assign(n, UNASSIGNED);
  std::size_t blocks = (num_vars + 63) / 64;
  adj_mat.assign(num_vars, std::vector<uint64_t>(blocks, 0ULL));

  // identity mapping
  for (std::size_t i = 0; i < num_vars; ++i) {
    g2l_map[i] = i;
    var_ids[i] = i;
  }
}

void PDAGwithAdjMat::set_as_complete() {
  const std::size_t blocks = (num_vars + 63) / 64;
  for (std::size_t i = 0; i < num_vars; ++i) {
    for (std::size_t b = 0; b < blocks; ++b) {
      uint64_t mask = (b == blocks - 1 && (num_vars % 64))
                          ? ((1ULL << (num_vars % 64)) - 1)
                          : ~0ULL;
      adj_mat[i][b] = mask;
    }
    _clr_arc(i, i);  // no self-loop
  }
}

PDAG PDAGwithAdjMat::to_pdag() const {
  PDAG p(num_vars);
  for (std::size_t vL = 0; vL < num_vars; ++vL) {
    for (auto uL : predecessorsL(vL)) {
      p.add_edge(var_ids[uL], var_ids[vL]);
      if (_has_arc(vL, uL)) p.add_edge(var_ids[vL], var_ids[uL]);
    }
  }
  return p;
}

PDAGwithAdjMat PDAGwithAdjMat::induced_subgraph(
    const PDAGwithAdjMat& src,
    const std::vector<size_t>& global_var_ids) {
  PDAGwithAdjMat res(global_var_ids.size());
  res.g2l_map.assign(src.num_vars, UNASSIGNED);
  res.var_ids.assign(res.num_vars, UNASSIGNED);

  for (std::size_t i = 0; i < res.num_vars; ++i) {
    res.g2l_map[global_var_ids[i]] = i;
    res.var_ids[i] = global_var_ids[i];
  }

  std::size_t blocks = (res.num_vars + 63) / 64;
  res.adj_mat.assign(res.num_vars, std::vector<uint64_t>(blocks, 0ULL));
  for (std::size_t i = 0; i < res.num_vars; ++i) {
    for (auto g_nb : src.undirected_neighborsL(global_var_ids[i])) {
      size_t j = res.g2l_map[g_nb];
      if (j == UNASSIGNED) continue;  // not in the subgraph
      res._set_arc(i, j);
      res._set_arc(j, i);
    }
  }
  return res;
}

bool PDAGwithAdjMat::_has_arc(std::size_t uL, std::size_t vL) const {
  _bounds(uL, num_vars);
  _bounds(vL, num_vars);
  std::size_t b = uL / 64, s = uL % 64;
  return (adj_mat[vL][b] & (1ULL << s)) != 0ULL;
}
void PDAGwithAdjMat::_set_arc(std::size_t uL, std::size_t vL) {
  _bounds(uL, num_vars);
  _bounds(vL, num_vars);
  std::size_t b = uL / 64, s = uL % 64;
  adj_mat[vL][b] |= (1ULL << s);
}
void PDAGwithAdjMat::_clr_arc(std::size_t uL, std::size_t vL) {
  _bounds(uL, num_vars);
  _bounds(vL, num_vars);
  std::size_t b = uL / 64, s = uL % 64;
  adj_mat[vL][b] &= ~(1ULL << s);
}

std::vector<std::size_t> PDAGwithAdjMat::predecessorsL(std::size_t vL) const {
  _bounds(vL, num_vars);
  std::vector<std::size_t> out;
  std::size_t blocks = (num_vars + 63) / 64;
  for (std::size_t b = 0; b < blocks; ++b) {
    uint64_t bits = adj_mat[vL][b];
    while (bits) {
      std::size_t s = __builtin_ctzll(bits);
      std::size_t uL = b * 64 + s;
      if (uL < num_vars) out.push_back(uL);
      bits &= bits - 1;
    }
  }
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::parentsL(std::size_t vL) const {
  std::vector<std::size_t> out;
  for (auto uL : predecessorsL(vL)) {
    if (!_has_arc(vL, uL)) out.push_back(uL);
  }
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighborsL(
    std::size_t vL) const {
  std::vector<std::size_t> out;
  for (auto uL : predecessorsL(vL)) {
    if (_has_arc(vL, uL)) out.push_back(uL);
  }
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors_withoutL(
    std::size_t vL,
    std::size_t excl) const {
  std::vector<std::size_t> out;
  for (auto uL : predecessorsL(vL)) {
    if (uL != excl && _has_arc(vL, uL)) out.push_back(uL);
  }
  return out;
}

/* TODO: why do we need g_all here? */
std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>>
PDAGwithAdjMat::decompose(const PDAGwithAdjMat& g_all) const {
  std::vector<size_t> desc_nodes_local;
  std::vector<size_t> desc_nodes;
  std::vector<size_t> remaining_nodes_local;
  for (std::size_t xL = 0; xL < num_vars; ++xL) {
    bool has_child = false;
    for (std::size_t yL = 0; yL < num_vars; ++yL) {
      if (has_directed_edgeL(xL, yL)) {
        has_child = true;
        break;
      }
    }
    if (has_child) {
      remaining_nodes_local.push_back(xL);
    } else {
      desc_nodes_local.push_back(xL);
      desc_nodes.push_back(var_ids[xL]);
    }
  }

  std::vector<std::vector<size_t>> asc_nodes_list;
  std::vector<char> seen(num_vars, 0);
  for (auto sL : remaining_nodes_local) {
    if (seen[sL]) continue;
    std::vector<size_t> asc_nodes;
    std::queue<size_t> q;
    q.push(sL);
    seen[sL] = 1;
    while (!q.empty()) {
      auto uL = q.front();
      q.pop();
      asc_nodes.push_back(var_ids[uL]);
      for (auto vL : remaining_nodes_local) {
        if (!seen[vL] && g_all.is_adjacentL(var_ids[uL], var_ids[vL])) {
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
// If (1) y-{x,z}, (2) x and z are non-adjacent, and (3) y is not in
// the sepset of x and z, then orient x->y<-z.
void PDAGwithAdjMat::orient_colliders(const Sepset& sepset) {
  for (std::size_t yL = 0; yL < num_vars; ++yL) {
    auto pres = predecessorsL(yL);  // {u | u->y or u<->y}
    for (std::size_t i = 0; i + 1 < pres.size(); ++i) {
      for (std::size_t j = i + 1; j < pres.size(); ++j) {
        auto xL = pres[i], zL = pres[j];
        if (is_adjacentL(xL, zL)) continue;
        size_t xG = var_ids[xL], yG = var_ids[yL], zG = var_ids[zL];
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

    /* ---------- Rule 1 ---------- */
    // If (1) x->y, (2) y-z, and (3) x and z are non-adjacent, then orient y->z.
    // Reason: orient y<-z would induce an unexpected unshielded collider
    // x->y<-z.
    for (std::size_t yL = 0; yL < num_vars; ++yL) {
      auto xL_candidates = parentsL(yL);
      auto zL_candidates = undirected_neighborsL(yL);
      for (auto xL : xL_candidates) {
        for (auto zL : zL_candidates) {
          if (!is_adjacentL(xL, zL)) {
            orient_edgeL(yL, zL);
            changed = true;
          }
        }
      }
    }

    /* ---------- Rule 2 ---------- */
    // If (1) x->y, (2) y->z, and (3) x-z, then orient x->z.
    // Reason: orient z<-x would induce an unexpected cycle x->y->z->x.
    for (std::size_t zL = 0; zL < num_vars; ++zL) {
      auto xL_candidates = undirected_neighborsL(zL);
      auto yL_candidates = parentsL(zL);
      for (auto xL : xL_candidates) {
        for (auto yL : yL_candidates) {
          if (has_directed_edgeL(xL, yL)) {
            orient_edgeL(xL, zL);
            changed = true;
          }
        }
      }
    }

    /* ---------- Rule 3 ---------- */
    // If (1) x-{y,z,w}, (2) {y,z} -> w, and (3) y and z are non-adjacent,
    // then orient x->w.
    // Reason: if we instead orient w->x, then to avoid cycles we would be
    // forced to orient {y,z}->x. But since y and z are non-adjacent, this would
    // create an unshielded collider y->x<-z.
    for (std::size_t wL = 0; wL < num_vars; ++wL) {
      auto xL_candidates = undirected_neighborsL(wL);
      auto yzL_candidates = parentsL(wL);
      for (std::size_t i = 0; i + 1 < yzL_candidates.size(); ++i) {
        for (std::size_t j = i + 1; j < yzL_candidates.size(); ++j) {
          auto yL = yzL_candidates[i], zL = yzL_candidates[j];
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
