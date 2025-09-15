#include "graph/pdag_with_adjmat.h"

#include <algorithm>
#include <limits>
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
    clr_arc(i, i);  // no self-loop
  }
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
    for (auto g_nb : src.undirected_neighbors(global_var_ids[i])) {
      size_t j = res.g2l_map[g_nb];
      if (j == UNASSIGNED) continue;  // not in the subgraph
      res.set_arc(i, j);
      res.set_arc(j, i);
    }
  }
  return res;
}

bool PDAGwithAdjMat::has_arc(std::size_t u, std::size_t v) const {
  _bounds(u, num_vars);
  _bounds(v, num_vars);
  std::size_t b = u / 64, s = u % 64;
  return (adj_mat[v][b] & (1ULL << s)) != 0ULL;
}
void PDAGwithAdjMat::set_arc(std::size_t u, std::size_t v) {
  _bounds(u, num_vars);
  _bounds(v, num_vars);
  std::size_t b = u / 64, s = u % 64;
  adj_mat[v][b] |= (1ULL << s);
}
void PDAGwithAdjMat::clr_arc(std::size_t u, std::size_t v) {
  _bounds(u, num_vars);
  _bounds(v, num_vars);
  std::size_t b = u / 64, s = u % 64;
  adj_mat[v][b] &= ~(1ULL << s);
}

std::vector<std::size_t> PDAGwithAdjMat::predecessors(std::size_t v) const {
  _bounds(v, num_vars);
  std::vector<std::size_t> out;
  std::size_t blocks = (num_vars + 63) / 64;
  for (std::size_t b = 0; b < blocks; ++b) {
    uint64_t bits = adj_mat[v][b];
    while (bits) {
      std::size_t s = __builtin_ctzll(bits);
      std::size_t u = b * 64 + s;
      if (u < num_vars) out.push_back(u);
      bits &= bits - 1;
    }
  }
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::parents(std::size_t v) const {
  std::vector<std::size_t> out;
  for (auto u : predecessors(v)) {
    if (!has_arc(v, u)) out.push_back(u);
  }
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors(
    std::size_t v) const {
  std::vector<std::size_t> out;
  for (auto u : predecessors(v)) {
    if (has_arc(v, u)) out.push_back(u);
  }
  return out;
}
std::vector<std::size_t> PDAGwithAdjMat::undirected_neighbors_without(
    std::size_t v,
    std::size_t excl) const {
  std::vector<std::size_t> out;
  for (auto u : predecessors(v)) {
    if (u != excl && has_arc(v, u)) out.push_back(u);
  }
  return out;
}

/* ---------- Rule V ---------- */
// If (1) y-{x,z}, (2) x and z are non-adjacent, and (3) y is not in
// the sepset of x and z, then orient x->y<-z.
void PDAGwithAdjMat::orient_colliders(const Sepset& sepset) {
  for (std::size_t y = 0; y < num_vars; ++y) {
    auto pres = predecessors(y);  // {u | u->y or u<->y}
    for (std::size_t i = 0; i + 1 < pres.size(); ++i) {
      for (std::size_t j = i + 1; j < pres.size(); ++j) {
        auto x = pres[i], z = pres[j];
        if (is_adjacent(x, z)) continue;
        size_t gx = var_ids[x], gy = var_ids[y], gz = var_ids[z];
        if (sepset[gx][gz].count(gy) == 0) {
          orient_edge(x, y);
          orient_edge(z, y);
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
    for (std::size_t y = 0; y < num_vars; ++y) {
      auto x_candidates = parents(y);
      auto z_candidates = undirected_neighbors(y);
      for (auto x : x_candidates) {
        for (auto z : z_candidates) {
          if (!is_adjacent(x, z)) {
            orient_edge(y, z);
            changed = true;
          }
        }
      }
    }

    /* ---------- Rule 2 ---------- */
    // If (1) x->y, (2) y->z, and (3) x-z, then orient x->z.
    // Reason: orient z<-x would induce an unexpected cycle x->y->z->x.
    for (std::size_t z = 0; z < num_vars; ++z) {
      auto x_candidates = undirected_neighbors(z);
      auto y_candidates = parents(z);
      for (auto x : x_candidates) {
        for (auto y : y_candidates) {
          if (has_directed_edge(x, y)) {
            orient_edge(x, z);
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
    for (std::size_t w = 0; w < num_vars; ++w) {
      auto x_candidates = undirected_neighbors(w);
      auto yz_candidates = parents(w);
      for (std::size_t i = 0; i + 1 < yz_candidates.size(); ++i) {
        for (std::size_t j = i + 1; j < yz_candidates.size(); ++j) {
          auto y = yz_candidates[i], z = yz_candidates[j];
          if (is_adjacent(y, z)) continue;
          for (auto x : x_candidates) {
            if (has_undirected_edge(x, y) && has_undirected_edge(x, z)) {
              orient_edge(x, w);
              changed = true;
            }
          }
        }
      }
    }
  }
}

std::vector<size_t> PDAGwithAdjMat::childless_nodes() const {
  std::vector<size_t> out;
  for (std::size_t x = 0; x < num_vars; ++x) {
    bool has_child = false;
    for (std::size_t y = 0; y < num_vars; ++y) {
      if (has_directed_edge(x, y)) {
        has_child = true;
        break;
      }
    }
    if (!has_child) out.push_back(var_ids[x]);
  }
  return out;
}

PDAG PDAGwithAdjMat::to_pdag() const {
  PDAG p(num_vars);
  for (std::size_t v = 0; v < num_vars; ++v) {
    for (auto u : predecessors(v)) {
      p.add_edge(var_ids[u], var_ids[v]);
      if (has_arc(v, u)) p.add_edge(var_ids[v], var_ids[u]);
    }
  }
  return p;
}