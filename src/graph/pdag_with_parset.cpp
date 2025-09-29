#include "graph/pdag_with_parset.h"

#include <algorithm>
#include <queue>
#include <stdexcept>

PDAGwithParSet::PDAGwithParSet(std::size_t n)
    : num_local_vars(n), num_global_vars(n) {
  gid2lid.assign(num_global_vars, UNASSIGNED);
  lid2gid.assign(num_local_vars, UNASSIGNED);
  in_parents.assign(num_local_vars, {});

  for (std::size_t i = 0; i < num_local_vars; ++i) {
    gid2lid[i] = i;
    lid2gid[i] = i;
  }
}

PDAGwithParSet PDAGwithParSet::induced_subgraph(const PDAGwithParSet& src,
                                                const std::vector<size_t>& S) {
  PDAGwithParSet res(S.size());
  res.num_global_vars = src.num_global_vars;
  res.gid2lid.assign(res.num_global_vars, UNASSIGNED);
  res.lid2gid.assign(res.num_local_vars, UNASSIGNED);
  for (std::size_t i = 0; i < res.num_local_vars; ++i) {
    const size_t g = S[i];
    res.gid2lid[g] = i;
    res.lid2gid[i] = g;
  }
  res.in_parents.assign(res.num_local_vars, {});

  for (std::size_t i = 0; i < res.num_local_vars; ++i) {
    const size_t g_i = res.lid2gid[i];
    for (auto g_nb : src.undirected_neighbors(g_i)) {
      const size_t j = res._to_local(g_nb);
      if (j == UNASSIGNED) continue;
      res.in_parents[i].insert(j);  // j -> i
      res.in_parents[j].insert(i);  // i -> j
    }
  }
  return res;
}

void PDAGwithParSet::set_as_complete() {
  for (std::size_t i = 0; i < num_local_vars; ++i) {
    for (std::size_t j = 0; j < num_local_vars; ++j) {
      if (i == j) continue;
      in_parents[i].insert(j);  // j -> i
    }
  }
}

std::unordered_set<std::size_t> PDAGwithParSet::predecessorsL(
    std::size_t vL) const {
  _bounds(vL, num_local_vars);
  return in_parents[vL];
}

std::unordered_set<std::size_t> PDAGwithParSet::parentsL(std::size_t vL) const {
  _bounds(vL, num_local_vars);
  std::unordered_set<size_t> out;
  for (auto uL : in_parents[vL]) {
    if (!in_parents[uL].count(vL)) out.insert(uL);
  }
  return out;
}

std::unordered_set<std::size_t> PDAGwithParSet::undirected_neighborsL(
    std::size_t vL) const {
  _bounds(vL, num_local_vars);
  std::unordered_set<size_t> out;
  for (auto uL : in_parents[vL]) {
    if (in_parents[uL].count(vL)) out.insert(uL);
  }
  return out;
}

std::unordered_set<std::size_t> PDAGwithParSet::undirected_neighbors_withoutL(
    std::size_t vL,
    std::size_t exclL) const {
  _bounds(vL, num_local_vars);
  std::unordered_set<size_t> out;
  for (auto uL : in_parents[vL]) {
    if (uL != exclL && in_parents[uL].count(vL)) out.insert(uL);
  }
  return out;
}

std::unordered_set<std::size_t> PDAGwithParSet::predecessors(
    std::size_t vG) const {
  const auto vL = _to_local(vG);
  if (vL == UNASSIGNED) throw std::runtime_error("v not in this graph");
  std::unordered_set<size_t> out;
  for (auto uL : predecessorsL(vL)) out.insert(lid2gid[uL]);
  return out;
}

std::unordered_set<std::size_t> PDAGwithParSet::parents(std::size_t vG) const {
  const auto vL = _to_local(vG);
  if (vL == UNASSIGNED) throw std::runtime_error("v not in this graph");
  std::unordered_set<size_t> out;
  for (auto uL : parentsL(vL)) out.insert(lid2gid[uL]);
  return out;
}

std::unordered_set<std::size_t> PDAGwithParSet::undirected_neighbors(
    std::size_t vG) const {
  const auto vL = _to_local(vG);
  if (vL == UNASSIGNED) throw std::runtime_error("v not in this graph");
  std::unordered_set<size_t> out;
  for (auto uL : undirected_neighborsL(vL)) out.insert(lid2gid[uL]);
  return out;
}

std::unordered_set<std::size_t> PDAGwithParSet::undirected_neighbors_without(
    std::size_t vG,
    std::size_t exclG) const {
  const auto vL = _to_local(vG);
  if (vL == UNASSIGNED) throw std::runtime_error("v not in this graph");
  const auto exclL = _to_local(exclG);
  std::unordered_set<size_t> out;
  for (auto uL : undirected_neighbors_withoutL(vL, exclL))
    out.insert(lid2gid[uL]);
  return out;
}

bool PDAGwithParSet::has_directed_edge(std::size_t uG, std::size_t vG) const {
  const auto uL = _to_local(uG), vL = _to_local(vG);
  if (uL == UNASSIGNED || vL == UNASSIGNED) return false;
  return has_directed_edgeL(uL, vL);
}

bool PDAGwithParSet::has_undirected_edge(std::size_t uG, std::size_t vG) const {
  const auto uL = _to_local(uG), vL = _to_local(vG);
  if (uL == UNASSIGNED || vL == UNASSIGNED) return false;
  return has_undirected_edgeL(uL, vL);
}

bool PDAGwithParSet::is_adjacent(std::size_t uG, std::size_t vG) const {
  const auto uL = _to_local(uG), vL = _to_local(vG);
  if (uL == UNASSIGNED || vL == UNASSIGNED) return false;
  return is_adjacentL(uL, vL);
}

void PDAGwithParSet::remove_undirected_edge(std::size_t uG, std::size_t vG) {
  const auto uL = _to_local(uG), vL = _to_local(vG);
  if (uL == UNASSIGNED || vL == UNASSIGNED) return;
  remove_undirected_edgeL(uL, vL);
}

void PDAGwithParSet::orient_edge(std::size_t uG, std::size_t vG) {
  const auto uL = _to_local(uG), vL = _to_local(vG);
  if (uL == UNASSIGNED || vL == UNASSIGNED) return;
  orient_edgeL(uL, vL);
}

PDAG PDAGwithParSet::to_pdag() const {
  PDAG p(num_global_vars);
  for (std::size_t vL = 0; vL < num_local_vars; ++vL) {
    const auto vG = lid2gid[vL];
    for (auto uL : predecessorsL(vL)) {
      const auto uG = lid2gid[uL];
      p.add_edge(uG, vG);
      if (in_parents[uL].count(vL)) p.add_edge(vG, uG);
    }
  }
  return p;
}

/* ---------- Rule V ----------
 * If (1) y-{x,z}, (2) x and z are non-adjacent, and
 * (3) y is not in the sepset of x and z, then orient x->y<-z.
 */
void PDAGwithParSet::orient_colliders(const Sepset& sepset) {
  for (std::size_t yL = 0; yL < num_local_vars; ++yL) {
    const auto& pres = in_parents[yL];  // {u | u->y or u<->y}
    for (auto it1 = pres.begin(); it1 != pres.end(); ++it1) {
      for (auto it2 = std::next(it1); it2 != pres.end(); ++it2) {
        const auto xL = *it1, zL = *it2;
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

void PDAGwithParSet::apply_meeks_rules() {
  bool changed = true;
  while (changed) {
    changed = false;

    /* ---------- Rule 1 ----------
     * If (1) x->y, (2) y-z, and (3) x and z are non-adjacent, then orient y->z.
     */
    for (std::size_t yL = 0; yL < num_local_vars; ++yL) {
      auto xLs = parentsL(yL);
      auto zLs = undirected_neighborsL(yL);
      for (auto xL : xLs) {
        for (auto zL : zLs) {
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
      auto xLs = undirected_neighborsL(zL);
      auto yLs = parentsL(zL);
      for (auto xL : xLs) {
        for (auto yL : yLs) {
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
      auto xLs = undirected_neighborsL(wL);
      auto yzLs = parentsL(wL);
      for (auto it1 = yzLs.begin(); it1 != yzLs.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != yzLs.end(); ++it2) {
          const auto yL = *it1, zL = *it2;
          if (is_adjacentL(yL, zL)) continue;
          for (auto xL : xLs) {
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

std::pair<std::vector<size_t>, std::vector<std::vector<size_t>>>
PDAGwithParSet::decompose(const PDAGwithParSet& g_all) const {
  std::vector<size_t> desc_nodes_local, desc_nodes, remaining_nodes_local;

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
