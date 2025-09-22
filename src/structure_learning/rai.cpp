#include "structure_learning/rai.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <numeric>
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
  std::vector<size_t> vars = cand_conds;
  vars.push_back(x);
  vars.push_back(y);
  std::sort(vars.begin(), vars.end());
  ContingencyTable<false> ct(vars, df);
  for (auto& Z : gen_combs(cand_conds, order_n)) {
    // std::vector<size_t> vars = Z;
    // vars.push_back(x);
    // vars.push_back(y);
    // std::sort(vars.begin(), vars.end());
    // ContingencyTable<false> ct(vars, df);
    if (citest<false>(x, y, Z, ct, test)) {
      g_all.remove_undirected_edgeL(x, y);
      if (g_sub) {
        int lx = g_sub->g2l_map[x], ly = g_sub->g2l_map[y];
        if (lx >= 0 && ly >= 0)
          g_sub->remove_undirected_edgeL((size_t)lx, (size_t)ly);
      }
      sepset[x][y].insert(Z.begin(), Z.end());
      sepset[y][x] = sepset[x][y];
      return true;
    }
  }
  return false;
}

static std::vector<size_t> extract_parents_from_exo(
    size_t v,
    const PDAGwithAdjMat& g,
    const std::vector<char>& in_exo) {
  std::vector<size_t> pa_exo;
  for (auto u : g.undirected_neighborsL(v))
    if (in_exo[u]) pa_exo.push_back(u);
  return pa_exo;
}

static std::vector<size_t> extract_potential_parents_from_sub(
    size_t v,
    const PDAGwithAdjMat& g_sub) {
  std::vector<size_t> pap_sub;
  for (auto uL : g_sub.predecessorsL(g_sub.g2l_map[v]))
    pap_sub.push_back(g_sub.var_ids[uL]);
  return pap_sub;
}

struct Subproblem {
  size_t order_n;
  std::vector<size_t> sub_nodes;
  std::vector<size_t> exo_nodes;
  PDAGwithAdjMat g_prev;
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
  const PDAGwithAdjMat& g_prev = curr.g_prev;

  PDAGwithAdjMat g_sub =
      PDAGwithAdjMat::induced_subgraph(g_all, curr.sub_nodes);
  auto in_exo = as_membership(curr.exo_nodes, g_all.num_vars);

  /* Exit condition */
  if (curr.order_n > ctx.max_cond_vars) return;
  bool has_sepset_candidates = false;
  for (auto y : curr.sub_nodes) {
    std::vector<size_t> pa_exo = extract_parents_from_exo(y, g_all, in_exo);
    std::vector<size_t> pap_sub = extract_potential_parents_from_sub(y, g_sub);
    auto base = merge_vectors(pa_exo, pap_sub);
    if (base.size() >= curr.order_n + 1) {
      has_sepset_candidates = true;
      break;
    }
  }
  if (!has_sepset_candidates) return;

  /* Stage A*/
  for (auto y : curr.sub_nodes) {
    std::vector<size_t> pa_exo = extract_parents_from_exo(y, g_all, in_exo);
    std::vector<size_t> pap_sub = extract_potential_parents_from_sub(y, g_sub);
    auto base = merge_vectors(pa_exo, pap_sub);

    for (auto x : pa_exo) {
      std::vector<size_t> cand;
      cand.reserve(base.size());
      for (auto v : base)
        if (v != x) cand.push_back(v);
      try_remove_and_apply(
          x, y, cand, curr.order_n, ctx.df, ctx.test, g_all, nullptr, sepset);
    }
  }
  g_sub.orient_colliders(sepset);
  g_sub.apply_meeks_rules();

  /* Stage B */
  for (auto y : curr.sub_nodes) {
    std::vector<size_t> pa_exo = extract_parents_from_exo(y, g_all, in_exo);
    std::vector<size_t> pap_sub = extract_potential_parents_from_sub(y, g_sub);
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
  auto [desc_nodes, asc_nodes_list] = g_sub.decompose(g_all);

  /* Summary of current subproblem */
  // DEBUG("RAI order_n=" << curr.order_n << ", sub_nodes:{"
  //                      << fmt_vec(curr.sub_nodes) << "}"
  //                      << ", exo_nodes:{" << fmt_vec(curr.exo_nodes) << "}");
  // for (auto& asc_nodes : asc_nodes_list)
  //   DEBUG("  ├─ asc_nodes={" << fmt_vec(asc_nodes) << "}");
  // DEBUG("  └─ desc_nodes={" << fmt_vec(desc_nodes) << "}");

  /* Stage C */
  for (auto& asc_nodes : asc_nodes_list) {
    Subproblem next{curr.order_n + 1, asc_nodes, curr.exo_nodes, g_sub};
    rai_recursive(next, g_all, sepset, ctx);
  }

  /* Stage D */
  std::vector<size_t> exo_nodes_for_desc = curr.exo_nodes;
  for (auto& asc_nodes : asc_nodes_list) {
    for (auto v : asc_nodes) {
      exo_nodes_for_desc.push_back(v);
    }
  }
  std::sort(exo_nodes_for_desc.begin(), exo_nodes_for_desc.end());
  Subproblem next{curr.order_n + 1, desc_nodes, exo_nodes_for_desc, g_sub};
  rai_recursive(next, g_all, sepset, ctx);
}

PDAG rai(const DataframeWrapper& df,
         const CITestType& test,
         size_t max_cond_vars) {
  const size_t n = df.num_of_vars;
  RAIContext ctx{df, test, max_cond_vars};

  PDAGwithAdjMat g_all(n);
  g_all.set_as_complete();
  Sepset sepset(n, std::vector<std::unordered_set<size_t>>(n));

  std::vector<size_t> all_nodes(n);
  for (size_t i = 0; i < n; ++i) all_nodes[i] = i;
  Subproblem init{0, all_nodes, {}, g_all};

  // GO!
  rai_recursive(init, g_all, sepset, ctx);

  g_all.orient_colliders(sepset);
  g_all.apply_meeks_rules();
  return g_all.to_pdag();
}