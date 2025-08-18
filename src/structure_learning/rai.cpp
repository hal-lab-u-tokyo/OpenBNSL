#include "structure_learning/rai.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <numeric>
#include <queue>
#include <unordered_set>
#include <vector>

#include "base/contingency_table.h"
#include "base/dataframe_wrapper.h"
#include "citest/citest.h"
#include "graph/pdag_with_adjmat.h"
#include "utils/gen_comb.h"

using Sepset = std::vector<std::vector<std::unordered_set<size_t>>>;

/* ============================================================
 *  Utilities
 * ============================================================ */

/* --- Collider orientation (PC 実装相当) --- */
static void orient_colliders(PDAGwithAdjMat& g, const Sepset& sepset) {
  const size_t n = g.num_vars;
  for (size_t y = 0; y < n; ++y) {
    bool changed = true;
    while (changed) {
      changed = false;
      auto neigh = g.undirected_neighbors(y);
      for (size_t i = 0; i + 1 < neigh.size(); ++i) {
        for (size_t j = i + 1; j < neigh.size(); ++j) {
          size_t x = neigh[i], z = neigh[j];
          if (g.is_adjacent(x, z)) continue;
          if (sepset[x][z].count(y) == 0) {
            g.orient_edge(x, y);
            g.orient_edge(z, y);
            changed = true;
          }
        }
      }
    }
  }
}

// /* --- bitset 的 membership ベクタ --- */
// static std::vector<char> as_membership(const std::vector<size_t>& nodes,
//                                        size_t nvars) {
//   std::vector<char> in(nvars, 0);
//   for (auto v : nodes) in[v] = 1;
//   return in;
// }

// /* --- 部分集合 in に含まれるものだけ抽出 --- */
// static std::vector<size_t> intersect_only(const std::vector<size_t>& xs,
//                                           const std::vector<char>& in) {
//   std::vector<size_t> out;
//   out.reserve(xs.size());
//   for (auto v : xs)
//     if (in[v]) out.push_back(v);
//   return out;
// }

// /* Pap(Y, Gstart) = Adj(Y) \ Ch(Y) を Gstart 制限で取得 */
// static std::vector<size_t> potential_parents_in_subgraph(
//     const PDAGwithAdjMat& g,
//     size_t y,
//     const std::vector<char>& in_start) {
//   std::vector<size_t> res;
//   for (auto u : g.neighbors(y)) {
//     if (!in_start[u]) continue;
//     bool is_child = g.has_directed_edge(y, u) && !g.has_directed_edge(u, y);
//     if (!is_child) res.push_back(u);  // Adj \ Ch
//   }
//   return res;
// }

// /* Pa(Y, Gex) = { X∈Gex | X → Y } */
// static std::vector<size_t> parents_from_exo(const PDAGwithAdjMat& g,
//                                             size_t y,
//                                             const std::vector<char>& in_exo)
//                                             {
//   std::vector<size_t> res;
//   for (size_t x = 0; x < g.num_vars; ++x) {
//     if (!in_exo[x]) continue;
//     if (g.has_directed_edge(x, y) && !g.has_directed_edge(y, x))
//       res.push_back(x);
//   }
//   return res;
// }

// /* Gstart 内で directed children を持たないノード（最下位トポロジの集合） */
// static std::vector<size_t> lowest_topology_nodes(
//     const PDAGwithAdjMat& g,
//     const std::vector<char>& in_start) {
//   std::vector<size_t> sinks;
//   for (size_t v = 0; v < g.num_vars; ++v) {
//     if (!in_start[v]) continue;
//     bool has_child = false;
//     for (auto c : g.directed_children(v)) {
//       if (in_start[c]) {
//         has_child = true;
//         break;
//       }
//     }
//     if (!has_child) sinks.push_back(v);
//   }
//   return sinks;
// }

// /* in_subset に制限した無向連結成分（向きは無視） */
// static std::vector<std::vector<size_t>> connected_components_in_subset(
//     const PDAGwithAdjMat& g,
//     const std::vector<char>& in_subset) {
//   std::vector<char> vis(g.num_vars, 0);
//   std::vector<std::vector<size_t>> comps;
//   for (size_t s = 0; s < g.num_vars; ++s) {
//     if (!in_subset[s] || vis[s]) continue;
//     std::vector<size_t> comp;
//     std::queue<size_t> q;
//     q.push(s);
//     vis[s] = 1;
//     while (!q.empty()) {
//       size_t u = q.front();
//       q.pop();
//       comp.push_back(u);
//       for (auto w : g.neighbors(u)) {
//         if (!in_subset[w] || vis[w]) continue;
//         vis[w] = 1;
//         q.push(w);
//       }
//     }
//     comps.push_back(std::move(comp));
//   }
//   return comps;
// }

// /* |S|=order_n の CI を試し、独立が見つかれば辺を削除＋sepset 記録（true
//  * で削除済） */
// static bool try_separate_and_remove_edge(size_t x,
//                                          size_t y,
//                                          const std::vector<size_t>&
//                                          cand_conds, size_t order_n, const
//                                          DataframeWrapper& df, const
//                                          CITestType& test, PDAGwithAdjMat& g,
//                                          Sepset& sepset) {
//   if (cand_conds.size() < order_n) return false;

//   for (auto& Z : gen_combs(cand_conds, order_n)) {
//     std::vector<size_t> vars = Z;
//     vars.push_back(x);
//     vars.push_back(y);
//     std::sort(vars.begin(), vars.end());

//     ContingencyTable<true> ct(vars, df);
//     if (citest<true>(x, y, Z, ct, test)) {
//       // 独立 → 辺を除去（向きに応じて削除）
//       if (g.has_undirected_edge(x, y)) {
//         g.remove_undirected_edge(x, y);
//       } else if (g.has_directed_edge(x, y)) {
//         g.remove_directed_edge(x, y);
//       } else if (g.has_directed_edge(y, x)) {
//         g.remove_directed_edge(y, x);
//       }
//       sepset[x][y].insert(Z.begin(), Z.end());
//       sepset[y][x] = sepset[x][y];
//       return true;
//     }
//   }
//   return false;
// }

// /* Exit 条件：∀Y∈Gstart, |Pa(Y,Gex) ∪ Pap(Y,Gstart)| < n+1 なら終了 */
// static bool exit_condition(const PDAGwithAdjMat& g,
//                            const std::vector<size_t>& Gstart,
//                            const std::vector<size_t>& Gex,
//                            size_t order_n) {
//   auto in_start = as_membership(Gstart, g.num_vars);
//   auto in_exo = as_membership(Gex, g.num_vars);

//   for (auto y : Gstart) {
//     auto pa_ex = parents_from_exo(g, y, in_exo);
//     auto pap_in = potential_parents_in_subgraph(g, y, in_start);

//     std::sort(pa_ex.begin(), pa_ex.end());
//     std::sort(pap_in.begin(), pap_in.end());

//     std::vector<size_t> uni;
//     std::set_union(pa_ex.begin(),
//                    pa_ex.end(),
//                    pap_in.begin(),
//                    pap_in.end(),
//                    std::back_inserter(uni));

//     if (uni.size() >= order_n + 1) return false;
//   }
//   return true;
// }

// /* ============================================================
//  *  RAI recursion (Stages A → B → C → D @ order n)
//  * ============================================================ */
// static void rai_recursive(size_t order_n,
//                           size_t max_cond_vars,
//                           const DataframeWrapper& df,
//                           const CITestType& test,
//                           bool apply_meek_r4,
//                           std::vector<size_t> Gstart,
//                           std::vector<size_t> Gex,
//                           PDAGwithAdjMat& gall,
//                           Sepset& sepset) {
//   if (order_n > max_cond_vars) return;
//   // if (Gstart.empty()) return;
//   // if (exit_condition(gall, Gstart, Gex, order_n)) return;

//   auto in_start = as_membership(Gstart, gall.num_vars);
//   auto in_exo = as_membership(Gex, gall.num_vars);

//   /* ---- Stage A: thin link G_{\text{start}} <- G_{\text{ex}} ---- */
//   // \forall Y \in G_{\text{start}}, X \in Pa(y, G_{\text{ex}})
//   // Remove edge (X, Y)
//   // if \exists S \subseteq \left{
//   //    Pa_\text{p}(Y,G_{\text{start}})
//   //    \cup Pa(Y,G_{\text{ex}})
//   //    \setminus X
//   // \right} s.t. |S| = n_\text{CI} and X \perp\!\!\!\perp Y | S
//   for (auto y : Gstart) {
//     auto pa_ex = parents_from_exo(gall, y, in_exo);
//     auto pap_in = potential_parents_in_subgraph(gall, y, in_start);
//     // if (pa_ex.empty()) continue;
//     std::vector<size_t> base = pa_ex;
//     base.insert(base.end(), pap_in.begin(), pap_in.end());
//     std::sort(base.begin(), base.end());
//     base.erase(std::unique(base.begin(), base.end()), base.end());

//     for (auto x : pa_ex) {
//       if (!gall.is_adjacent(x, y)) continue;
//       std::vector<size_t> cand;
//       cand.reserve(base.size());
//       for (auto v : base)
//         if (v != x) cand.push_back(v);
//       try_separate_and_remove_edge(x, y, cand, order_n, df, test, gall,
//       sepset);
//     }
//   }
//   orient_colliders(gall, sepset);
//   // gall.apply_meeks_rules(apply_meek_r4);
//   /* ---- Stage A ---- */

//   /* ---- Stage B: thin/orient/decompose G_{\text{start}} ---- */
//   // \forall Y \in G_{\text{start}}, X \in Pa(Y, G_{\text{start}})
//   // Remove edge (X, Y)
//   // if \exists S \subseteq \left{
//   //    Pa_\text{p}(Y,G_{\text{start}})
//   //    \cup Pa(Y,G_{\text{ex}})
//   //    \setminus X
//   // \right} s.t. |S| = n_\text{CI} and X \perp\!\!\!\perp Y | S

//   for (auto y : Gstart) {
//     auto pa_ex = parents_from_exo(gall, y, in_exo);
//     auto pap_in = potential_parents_in_subgraph(gall, y, in_start);
//     // if (pap_in.empty()) continue;

//     std::vector<size_t> base = pa_ex;
//     base.insert(base.end(), pap_in.begin(), pap_in.end());
//     std::sort(base.begin(), base.end());
//     base.erase(std::unique(base.begin(), base.end()), base.end());

//     for (auto x : pap_in) {
//       if (!gall.is_adjacent(x, y)) continue;
//       std::vector<size_t> cand;
//       cand.reserve(base.size());
//       for (auto v : base)
//         if (v != x) cand.push_back(v);
//       try_separate_and_remove_edge(x, y, cand, order_n, df, test, gall,
//       sepset);
//     }
//   }

//   orient_colliders(gall, sepset);
//   // gall.apply_meeks_rules(apply_meek_r4);

//   // B3: 最下位トポロジのノード群を GD（descendants）に
//   auto GD = lowest_topology_nodes(gall, in_start);
//   if (GD.empty()) {
//     // フォールバック：最小サイズの連結成分を GD に
//     auto comps = connected_components_in_subset(gall, in_start);
//     if (!comps.empty()) {
//       size_t best_i = 0, best_sz = comps[0].size();
//       for (size_t i = 1; i < comps.size(); ++i) {
//         if (comps[i].size() < best_sz) {
//           best_sz = comps[i].size();
//           best_i = i;
//         }
//       }
//       GD = comps[best_i];
//     }
//   }

//   // B4: GA1..GAk（ancestors）＝ Gstart \ GD の連結成分
//   auto in_gd = as_membership(GD, gall.num_vars);
//   std::vector<char> in_anc = in_start;
//   for (auto v : GD) in_anc[v] = 0;
//   auto GAs = connected_components_in_subset(gall, in_anc);

//   // GA → GD の未向き辺を明示的に向ける（トポロジ順）
//   for (auto& GA : GAs) {
//     auto in_ga = as_membership(GA, gall.num_vars);
//     for (size_t a = 0; a < gall.num_vars; ++a) {
//       if (!in_ga[a]) continue;
//       for (auto d : GD) {
//         if (gall.has_undirected_edge(a, d)) gall.orient_edge(a, d);
//       }
//     }
//   }
//   // gall.apply_meeks_rules(apply_meek_r4);
//   /* ---- Stage B ---- */

//   /* ---- Stage C: recurse on each GAi with n+1 ---- */
//   for (auto& GA : GAs) {
//     rai_recursive(order_n + 1,
//                   max_cond_vars,
//                   df,
//                   test,
//                   apply_meek_r4,
//                   GA,  /* Gstart */
//                   Gex, /* Gex unchanged */
//                   gall,
//                   sepset);
//   }
//   /* ---- Stage C ---- */

//   /* ---- Stage D: recurse on GD with n+1 (exogenous = GA∪Gex) ---- */
//   std::vector<size_t> GexD = Gex;
//   for (auto& GA : GAs) {
//     GexD.insert(GexD.end(), GA.begin(), GA.end());
//   }
//   std::sort(GexD.begin(), GexD.end());
//   GexD.erase(std::unique(GexD.begin(), GexD.end()), GexD.end());

//   rai_recursive(order_n + 1,
//                 max_cond_vars,
//                 df,
//                 test,
//                 apply_meek_r4,
//                 GD,   /* Gstart = GD */
//                 GexD, /* Gex = GA1..GAk ∪ Gex */
//                 gall,
//                 sepset);
//   /* ---- Stage D ---- */
// }

/* ============================================================
 *  Public API
 * ============================================================ */
PDAG rai(const DataframeWrapper& df,
         const CITestType& test,
         size_t max_cond_vars,
         bool apply_meek_r4) {
  const size_t n = df.num_of_vars;
  PDAGwithAdjMat gall(n);
  gall.complete_graph();
  Sepset sepset(n, std::vector<std::unordered_set<size_t>>(n));

  // std::vector<size_t> all_nodes(n);
  // std::iota(all_nodes.begin(), all_nodes.end(), 0);
  // rai_recursive(/*order_n=*/0,
  //               max_cond_vars,
  //               df,
  //               test,
  //               apply_meek_r4,
  //               /*Gstart=*/all_nodes,
  //               /*Gex=*/{},
  //               gall,
  //               sepset);
  // orient_colliders(gall, sepset);
  // gall.apply_meeks_rules(apply_meek_r4);
  return gall.to_pdag();
}
