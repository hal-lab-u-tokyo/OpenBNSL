#include "structure_learning/rai.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

#include "base/contingency_table.h"
#include "base/dataframe_wrapper.h"
#include "citest/citest.h"
#include "graph/pdag_with_adjmat.h"
#include "utils/format.h"
#include "utils/gen_comb.h"
#include "utils/logging.h"

/* ============================================================
 *  Utilities
 * ============================================================ */
template <typename T>
static std::vector<T> merge_vectors(const std::vector<T>& a,
                                    const std::vector<T>& b) {
  std::vector<T> r;
  r.reserve(a.size() + b.size());
  r.insert(r.end(), a.begin(), a.end());
  r.insert(r.end(), b.begin(), b.end());
  std::sort(r.begin(), r.end());
  r.erase(std::unique(r.begin(), r.end()), r.end());
  return r;
}

static std::vector<char> as_membership(const std::vector<size_t>& nodes,
                                       size_t n) {
  std::vector<char> in(n, 0);
  for (auto v : nodes) in[v] = 1;
  return in;
}

/*
 * In phase A, we try to remove edges from only g_all
 * In phase B, we try to remove edges from g_all and g_sub
 */
static bool try_remove_and_apply(size_t x,
                                 size_t y,
                                 const std::vector<size_t>& cand_conds,
                                 size_t order_n,
                                 const DataframeWrapper& df,
                                 const CITestType& test,
                                 PDAGwithAdjMat& g_all,
                                 PDAGwithAdjMat* g_sub,
                                 Sepset& sepset) {
  if (cand_conds.size() < order_n) return false;
  if (!g_all.is_adjacent(x, y))
    throw std::logic_error("try_remove_and_apply for non-adjacent nodes");

  for (auto& Z : gen_combs(cand_conds, order_n)) {
    std::vector<size_t> vars = Z;
    vars.push_back(x);
    vars.push_back(y);
    std::sort(vars.begin(), vars.end());
    ContingencyTable<true> ct(vars, df);
    if (citest<true>(x, y, Z, ct, test)) {
      g_all.remove_undirected_edge(x, y);
      if (g_sub) {
        int lx = g_sub->g2l_map[x], ly = g_sub->g2l_map[y];
        if (lx >= 0 && ly >= 0)
          g_sub->remove_undirected_edge((size_t)lx, (size_t)ly);
      }
      sepset[x][y].insert(Z.begin(), Z.end());
      sepset[y][x] = sepset[x][y];
      return true;
    }
  }
  return false;
}

struct Subproblem {
  size_t order_n;
  std::vector<size_t> sub_nodes;
  std::vector<size_t> exo_nodes;
};

struct RAIContext {
  const DataframeWrapper& df;
  const CITestType& test;
  size_t max_cond_vars;
};

static void rai_recursive(const Subproblem& curr,
                          PDAGwithAdjMat& g_all,
                          Sepset& sepset,
                          const RAIContext& ctx) {
  if (curr.order_n > ctx.max_cond_vars) return;
  // TODO: check exit condition

  PDAGwithAdjMat g_sub =
      PDAGwithAdjMat::induced_subgraph(g_all, curr.sub_nodes);
  auto in_exo = as_membership(curr.exo_nodes, g_all.num_vars);

  /* Stage A*/
  for (auto y : curr.sub_nodes) {
    std::vector<size_t> pa_exo;
    for (auto u : g_all.undirected_neighbors(y))
      if (in_exo[u]) pa_exo.push_back(u);
    // for (auto u : g_sub.undirected_neighbors(g_sub.g2l_map[y]))
    //   if (in_exo[g_sub.var_ids[u]]) pa_exo.push_back(g_sub.var_ids[u]);

    std::vector<size_t> pap_sub;
    for (auto uL : g_sub.predecessors(g_sub.g2l_map[y]))
      pap_sub.push_back(g_sub.var_ids[uL]);

    auto base = merge_vectors(pa_exo, pap_sub);
    for (auto x : pa_exo) {
      std::vector<size_t> cand;
      cand.reserve(base.size());
      for (auto v : base)
        if (v != x) cand.push_back(v);
      try_remove_and_apply(x,
                           y,
                           cand,
                           curr.order_n,
                           ctx.df,
                           ctx.test,
                           g_all,
                           /*g_sub=*/nullptr,
                           sepset);
    }
  }
  g_sub.orient_colliders(sepset);
  g_sub.apply_meeks_rules();

  /* Stage B */
  for (auto y : curr.sub_nodes) {
    std::vector<size_t> pa_exo;
    for (auto u : g_all.undirected_neighbors(y))
      if (in_exo[u]) pa_exo.push_back(u);

    std::vector<size_t> pap_sub;
    for (auto uL : g_sub.predecessors(g_sub.g2l_map[y]))
      pap_sub.push_back(g_sub.var_ids[uL]);

    auto base = merge_vectors(pa_exo, pap_sub);
    for (auto x : pap_sub) {
      std::vector<size_t> cand;
      cand.reserve(base.size());
      for (auto v : base)
        if (v != x) cand.push_back(v);

      try_remove_and_apply(
          x, y, cand, curr.order_n, ctx.df, ctx.test, g_all, &g_sub, sepset);
    }
  }
  g_sub.orient_colliders(sepset);
  g_sub.apply_meeks_rules();

  /* Decomposition */
  std::vector<size_t> desc_nodes = g_sub.childless_nodes();
  std::vector<char> in_desc_nodes = as_membership(desc_nodes, g_all.num_vars);
  std::vector<size_t> remaining_nodes;
  for (auto v : curr.sub_nodes)
    if (!in_desc_nodes[v]) remaining_nodes.push_back(v);

  std::vector<std::vector<size_t>> asc_nodes_list;
  if (!remaining_nodes.empty()) {
    std::vector<char> seen(g_all.num_vars, 0);
    for (auto s : remaining_nodes) {
      if (seen[s]) continue;
      std::vector<size_t> comp;
      std::queue<size_t> q;
      q.push(s);
      seen[s] = 1;
      while (!q.empty()) {
        auto u = q.front();
        q.pop();
        comp.push_back(u);
        for (auto v : remaining_nodes) {
          if (!seen[v] && g_all.is_adjacent(u, v)) {
            seen[v] = 1;
            q.push(v);
          }
        }
      }
      asc_nodes_list.push_back(std::move(comp));
    }
  }

  std::vector<size_t> exo_nodes_for_desc = curr.exo_nodes;
  for (auto& asc_nodes : asc_nodes_list) {
    for (auto v : asc_nodes) {
      exo_nodes_for_desc.push_back(v);
    }
  }
  std::sort(exo_nodes_for_desc.begin(), exo_nodes_for_desc.end());
  exo_nodes_for_desc.erase(
      std::unique(exo_nodes_for_desc.begin(), exo_nodes_for_desc.end()),
      exo_nodes_for_desc.end());

  /* Summary of current subproblem */
  DEBUG("RAI order_n=" << curr.order_n << ", sub_nodes:{"
                       << fmt_vec(curr.sub_nodes) << "}"
                       << ", exo_nodes:{" << fmt_vec(curr.exo_nodes) << "}");
  for (auto& asc_nodes : asc_nodes_list)
    DEBUG("  ├─ asc_nodes={" << fmt_vec(asc_nodes) << "}");
  DEBUG("  └─ desc_nodes={" << fmt_vec(desc_nodes) << "}");

  /* Stage C */
  for (auto& asc_nodes : asc_nodes_list) {
    Subproblem sb_next{curr.order_n + 1, asc_nodes, curr.exo_nodes};
    rai_recursive(sb_next, g_all, sepset, ctx);
  }

  /* Stage D */
  Subproblem sb_next{curr.order_n + 1, desc_nodes, exo_nodes_for_desc};
  rai_recursive(sb_next, g_all, sepset, ctx);
}

PDAG rai(const DataframeWrapper& df,
         const CITestType& test,
         size_t max_cond_vars) {
  const size_t n = df.num_of_vars;
  RAIContext ctx{df, test, max_cond_vars};

  PDAGwithAdjMat g_all(n);
  g_all.set_as_complete();
  Sepset sepset(n, std::vector<std::unordered_set<size_t>>(n));

  std::vector<size_t> all(n);
  for (size_t i = 0; i < n; ++i) all[i] = i;
  Subproblem sp_init{0, all, {}};

  // GO!
  rai_recursive(sp_init, g_all, sepset, ctx);

  g_all.orient_colliders(sepset);
  g_all.apply_meeks_rules();
  return g_all.to_pdag();
}