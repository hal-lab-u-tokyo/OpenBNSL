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
#include "utils/merge_set.h"

static std::vector<char> as_membership(const std::vector<size_t>& nodes,
                                       size_t n) {
  std::vector<char> in(n, 0);
  for (auto v : nodes)
    if (v < n) in[v] = 1;
  return in;
}

/*
 * In phase A, we try to remove edges from only g_all
 * In phase B, we try to remove edges from g_all and g_sub
 */
static bool try_remove_and_apply(size_t xG,
                                 size_t yG,
                                 const std::vector<size_t>& cand_condsG,
                                 size_t order_n,
                                 const DataframeWrapper& df,
                                 const CITestType& test,
                                 PDAGwithAdjMat* g_sub,
                                 PDAGwithAdjMat& g_all,
                                 Sepset& sepset) {
  // B: Create contingency table for {xG,yG} U pap(yG) once
  // std::vector<size_t> vars = cand_condsG;
  // vars.push_back(xG);
  // vars.push_back(yG);
  // std::sort(vars.begin(), vars.end());
  // ContingencyTable ct(vars, df);

  for (const auto& Z : gen_combs(cand_condsG, order_n)) {
    // A: Create contingency table for {xG,yG} U Z for each Z
    std::vector<size_t> vars = Z;
    vars.push_back(xG);
    vars.push_back(yG);
    std::sort(vars.begin(), vars.end());
    ContingencyTable ct(vars, df);

    if (citest(xG, yG, Z, ct, test)) {
      g_all.remove_undirected_edge(xG, yG);
      if (g_sub) g_sub->remove_undirected_edge(xG, yG);
      sepset[xG][yG].insert(Z.begin(), Z.end());
      sepset[yG][xG] = sepset[xG][yG];
      return true;
    }
  }
  return false;
}

static std::vector<size_t> extract_parents_from_exo(
    size_t vG,
    const PDAGwithAdjMat& g_all,
    const std::vector<char>& in_exo) {
  std::vector<size_t> pa_exo;
  for (auto uG : g_all.undirected_neighbors(vG)) {
    if (uG < in_exo.size() && in_exo[uG]) pa_exo.push_back(uG);
  }
  return pa_exo;
}

static std::vector<size_t> extract_potential_parents_from_sub(
    size_t vG,
    const PDAGwithAdjMat& g_sub) {
  return g_sub.predecessors(vG);
}

struct Subproblem {
  size_t order_n;
  std::vector<size_t> sub_nodes;  // global ids
  std::vector<size_t> exo_nodes;  // global ids
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
  PDAGwithAdjMat g_sub =
      PDAGwithAdjMat::induced_subgraph(g_all, curr.sub_nodes);
  auto in_exo = as_membership(curr.exo_nodes, g_all.num_global_vars);

  /* Exit condition */
  if (curr.order_n > ctx.max_cond_vars) return;
  bool has_sepset_candidates = false;
  for (auto yG : curr.sub_nodes) {
    const auto pa_exo = extract_parents_from_exo(yG, g_all, in_exo);
    const auto pap_sub = extract_potential_parents_from_sub(yG, g_sub);
    const auto base = merge_set(pa_exo, pap_sub);
    if (base.size() >= curr.order_n + 1) {
      has_sepset_candidates = true;
      break;
    }
  }
  if (!has_sepset_candidates) return;

  /* Stage A */
  for (auto yG : curr.sub_nodes) {
    const auto pa_exo = extract_parents_from_exo(yG, g_all, in_exo);
    const auto pap_sub = extract_potential_parents_from_sub(yG, g_sub);
    const auto base = merge_set(pa_exo, pap_sub);

    for (auto xG : pa_exo) {
      std::vector<size_t> cand;
      cand.reserve(base.size());
      for (auto vG : base)
        if (vG != xG) cand.push_back(vG);

      try_remove_and_apply(
          xG, yG, cand, curr.order_n, ctx.df, ctx.test, nullptr, g_all, sepset);
    }
  }
  g_sub.orient_colliders(sepset);
  g_sub.apply_meeks_rules();

  /* Stage B */
  for (auto yG : curr.sub_nodes) {
    const auto pa_exo = extract_parents_from_exo(yG, g_all, in_exo);
    const auto pap_sub = extract_potential_parents_from_sub(yG, g_sub);
    const auto base = merge_set(pa_exo, pap_sub);

    for (auto xG : pap_sub) {
      std::vector<size_t> cand;
      cand.reserve(base.size());
      for (auto vG : base)
        if (vG != xG) cand.push_back(vG);

      try_remove_and_apply(
          xG, yG, cand, curr.order_n, ctx.df, ctx.test, &g_sub, g_all, sepset);
    }
  }
  g_sub.orient_colliders(sepset);
  g_sub.apply_meeks_rules();

  /* Decomposition */
  auto [desc_nodes, asc_nodes_list] = g_sub.decompose(g_all);

  /* Stage C */
  for (auto& asc_nodes : asc_nodes_list) {
    Subproblem next{curr.order_n + 1, asc_nodes, curr.exo_nodes};
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
  Subproblem next{curr.order_n + 1, desc_nodes, exo_nodes_for_desc};
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
  Subproblem init{0, all_nodes, {}};

  // GO!
  rai_recursive(init, g_all, sepset, ctx);

  g_all.orient_colliders(sepset);
  g_all.apply_meeks_rules();
  return g_all.to_pdag();
}
