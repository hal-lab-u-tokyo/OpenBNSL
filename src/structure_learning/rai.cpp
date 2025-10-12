#include "structure_learning/rai.h"

#include "base/contingency_table.h"
#include "base/dataframe_wrapper.h"
#include "citest/citest.h"
#include "graph/pdag_with_adjmat.h"
#include "utils/gen_comb.h"
#include "utils/logging.h"

#include <deque>

std::vector<size_t> merge_sets(
    const std::vector<size_t>& a,
    const std::vector<size_t>& b) {
  std::vector<size_t> res;
  res.reserve(a.size() + b.size());
  res.insert(res.end(), a.begin(), a.end());
  res.insert(res.end(), b.begin(), b.end());
  std::sort(res.begin(), res.end());
  // No need, but just in case
  res.erase(std::unique(res.begin(), res.end()), res.end());
  return res;
}

struct Update {
  size_t x, y;
  std::vector<size_t> Z;
};

static inline bool bit_test(const std::vector<uint64_t>& bits, size_t i) {
  return (bits[i / 64] >> (i % 64)) & 1ULL;
}
static inline void bit_set(std::vector<uint64_t>& bits, size_t i) {
  bits[i / 64] |= (1ULL << (i % 64));
}
static inline void bit_clear(std::vector<uint64_t>& bits, size_t i) {
  bits[i / 64] &= ~(1ULL << (i % 64));
}

static std::pair<std::vector<uint64_t>, std::vector<std::vector<uint64_t>>>
decompose(const PDAGwithAdjMat& g,
          const std::vector<uint64_t>& sub_bits,
          const std::vector<uint64_t>& /*exo_bits not used*/) {
  const size_t n = g.num_vars;
  const size_t blocks = (n + 63) / 64;
  const auto sub_nodes = g._extract_indices_from_bits(sub_bits);

  std::vector<char> is_parent_of_someone_in_sub(n, 0);
  for (auto v : sub_nodes) {
    for (auto p : g.parents_in(v, sub_bits)) is_parent_of_someone_in_sub[p] = 1;
  }

  // Only for debugging, set all bits to 1 in desc_bits
  // std::vector<uint64_t> _desc_bits(blocks, 0ULL);
  // for (auto v : sub_nodes) bit_set(_desc_bits, v);
  // std::vector<std::vector<uint64_t>> _asc_bits_list;
  // return {std::move(_desc_bits), std::move(_asc_bits_list)};

  // A node that has no child nodes in the sub and cannot reach any node having child nodes in the sub via undirected edges.
  // As a result, any edge from a desc to any asc (that is from sub to exo) cannot exist in any subsequent orientation.
  std::vector<uint64_t> desc_bits(blocks, 0ULL);
  std::vector<char> visited(n, 0);
  for (auto seed : sub_nodes) {
    if (is_parent_of_someone_in_sub[seed] || visited[seed]) continue;

    std::deque<size_t> q;
    std::vector<size_t> comp;
    q.push_back(seed);
    visited[seed] = 1;
    bool touches_parent = false;

    while (!q.empty()) {
      auto u = q.front(); q.pop_front();
      comp.push_back(u);
      if (is_parent_of_someone_in_sub[u]) touches_parent = true;

      for (auto nb : g.undirected_neighbors_in(u, sub_bits)) {
        if (!visited[nb]) {
          visited[nb] = 1;
          q.push_back(nb);
        }
      }
    }

    if (!touches_parent) {
      for (auto v : comp) {
        bit_set(desc_bits, v);
      }
    }
  }

  std::vector<std::vector<uint64_t>> asc_bits_list;
  std::vector<char> seen(n, 0);

  for (auto v : sub_nodes) {
    if (bit_test(desc_bits, v) || seen[v]) continue;

    std::vector<uint64_t> bits(blocks, 0ULL);
    std::deque<size_t> q;
    q.push_back(v);
    seen[v] = 1;
    bit_set(bits, v);

    while (!q.empty()) {
      auto u = q.front(); q.pop_front();

      for (size_t w = 0; w < n; ++w) {
        if (!bit_test(sub_bits, w)) continue;
        if (bit_test(desc_bits, w)) continue;
        if (seen[w]) continue;
        if (g.is_adjacent(u, w)) {
          seen[w] = 1;
          bit_set(bits, w);
          q.push_back(w);
        }
      }
    }

    asc_bits_list.push_back(std::move(bits));
  }

  return {std::move(desc_bits), std::move(asc_bits_list)};
}


static void rai_recursive(const size_t k,
                          std::vector<uint64_t> sub_bits,
                          std::vector<uint64_t> exo_bits,
                          PDAGwithAdjMat& g,
                          Sepset& sepset,
                          const DataframeWrapper& df,
                          const CITestType& test,
                          size_t max_cond_vars) {
  auto start_collect_neighs = std::chrono::high_resolution_clock::now();
  const size_t n = g.num_vars;
  const auto sub_nodes = g._extract_indices_from_bits(sub_bits);
  const auto exo_nodes = g._extract_indices_from_bits(exo_bits);
  std::vector<std::vector<size_t>> pa_exo(n);
  std::vector<std::vector<size_t>> pap_sub(n);
  std::vector<std::vector<size_t>> pap_all(n);
  for (auto y : sub_nodes) {
    // pa_exo[y] = g.parents_in(y, exo_bits);
    pa_exo[y] = g.potential_parents_in(y, exo_bits);
    pap_sub[y] = g.potential_parents_in(y, sub_bits);
    pap_all[y] = merge_sets(pa_exo[y], pap_sub[y]);
  }

  /* Exit condition */
  if (k > max_cond_vars) return;
  bool has_sepset_candidates = false;
  for (auto y : sub_nodes) {
    if (pap_all[y].size() >= k + 1) {
      has_sepset_candidates = true;
      break;
    }
  }
  if (!has_sepset_candidates) return;
  auto end_collect_neighs = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_collect_neighs = end_collect_neighs - start_collect_neighs;
  INFO("[RAI] order=" << k
        << ", time collecting neighbors=" << elapsed_collect_neighs.count() << "s");

  /* Stage A */
  auto start_stageA = std::chrono::high_resolution_clock::now();
  std::vector<std::pair<size_t,size_t>> stageA_pairs; // [(child, parent)]
  for (auto y : sub_nodes) {
    for (auto x : pa_exo[y]) {
      stageA_pairs.emplace_back(y, x);
    }
  }

  size_t stageA_num_citests = 0;
  std::vector<Update> stageA_updates;

#pragma omp parallel default(none) \
  shared(k, df, test, stageA_pairs, stageA_num_citests, stageA_updates, pap_all)
  {
    size_t local_num_citests = 0;
    std::vector<Update> local_updates;

#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < stageA_pairs.size(); ++i) {
      auto [y, x] = stageA_pairs[i];
      const auto& pap_of_y = pap_all[y];
      if (pap_of_y.size() < k + 1) continue;
      std::vector<size_t> pap_of_y_without_x;
      pap_of_y_without_x.reserve(pap_of_y.size() - 1);
      for (auto v : pap_of_y)
        if (v != x) pap_of_y_without_x.push_back(v);
      std::sort(pap_of_y_without_x.begin(), pap_of_y_without_x.end());

      for (const auto& Z : gen_combs(pap_of_y_without_x, k)) {
        std::vector<size_t> vars = Z;
        vars.push_back(x); 
        vars.push_back(y);
        std::sort(vars.begin(), vars.end());
        ContingencyTable ct(vars, df);
        local_num_citests++;

        if (citest(x, y, Z, ct, test)) {
          Update update{x, y, std::vector<size_t>(Z.begin(), Z.end())};
          local_updates.push_back(update);
          break;
        }
      }
    } // end omp for
#pragma omp critical
    {
      stageA_num_citests += local_num_citests;
      stageA_updates.insert(stageA_updates.end(), local_updates.begin(), local_updates.end());
    }
  } // end omp parallel

  for (const auto& u : stageA_updates) {
    g.remove_edge(u.x, u.y);
    sepset[u.x][u.y].insert(u.Z.begin(), u.Z.end());
    sepset[u.y][u.x] = sepset[u.x][u.y];
  }
  g.orient_colliders(sepset);
  g.apply_meeks_rules();

  auto end_stageA = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_a = end_stageA - start_stageA;
  INFO("[RAI] order=" << k
        << ", edges removed in stage A=" << stageA_updates.size()
        << "/" << stageA_num_citests
        << ", time=" << elapsed_a.count() << "s");

  /* Stage B */
  auto start_stageB = std::chrono::high_resolution_clock::now();
  for (auto y : sub_nodes) {
    // pa_exo[y] = g.parents_in(y, exo_bits);
    pa_exo[y] = g.potential_parents_in(y, exo_bits);
    pap_sub[y] = g.potential_parents_in(y, sub_bits);
    pap_all[y] = merge_sets(pa_exo[y], pap_sub[y]);
  }

  std::vector<std::pair<size_t,size_t>> stageB_pairs;
  std::unordered_set<uint64_t> seen;
  for (auto y : sub_nodes) {
    for (auto x : pap_sub[y]) {
      auto u = std::min(x, y); auto v = std::max(x, y);
      uint64_t key = (uint64_t(u) << 32) | uint64_t(v);
      if (seen.insert(key).second) stageB_pairs.emplace_back(u, v);
    }
  }

  size_t stageB_num_citests = 0;
  std::vector<Update> stageB_updates;
#pragma omp parallel default(none) \
  shared(k, df, test, stageB_pairs, stageB_num_citests, stageB_updates, pap_all)
  {
    size_t local_num_citests = 0;
    std::vector<Update> local_updates;
#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < stageB_pairs.size(); ++i) {
      auto [x, y] = stageB_pairs[i];
      // iterate on smaller one first
      if (pap_all[x].size() > pap_all[y].size()) std::swap(x, y); 
      for (auto [u, v] : std::array{std::pair{x, y}, std::pair{y, x}}) {
        const auto& pap_of_u = pap_all[u];
        if (pap_of_u.size() < k + 1) continue;
        std::vector<size_t> pap_of_u_without_v;
        pap_of_u_without_v.reserve(pap_of_u.size() - 1);
        for (auto w : pap_of_u)
          if (w != v) pap_of_u_without_v.push_back(w);
        std::sort(pap_of_u_without_v.begin(), pap_of_u_without_v.end());

        for (const auto& Z : gen_combs(pap_of_u_without_v, k)) {
          std::vector<size_t> vars = Z;
          vars.push_back(u);
          vars.push_back(v);
          std::sort(vars.begin(), vars.end());
          ContingencyTable ct(vars, df);
          local_num_citests++;

          if (citest(u, v, Z, ct, test)) {
            Update update{x, y, std::vector<size_t>(Z.begin(), Z.end())};
            local_updates.push_back(update);
            goto NEXT_PAIR;
          }
        }
      }
      NEXT_PAIR:;
    } // end omp for
#pragma omp critical
    {
      stageB_num_citests += local_num_citests;
      stageB_updates.insert(stageB_updates.end(), local_updates.begin(), local_updates.end());
    }
  } // end omp parallel
  
  for (const auto& u : stageB_updates) {
    g.remove_edge(u.x, u.y);
    sepset[u.x][u.y].insert(u.Z.begin(), u.Z.end());
    sepset[u.y][u.x] = sepset[u.x][u.y];
  }
  g.orient_colliders(sepset);
  g.apply_meeks_rules();

  auto end_stageB = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_b = end_stageB - start_stageB;
  INFO("[RAI] order=" << k
        << ", edges removed in stage B=" << stageB_updates.size()
        << "/" << stageB_num_citests
        << ", time=" << elapsed_b.count() << "s");
  
  /* Decomposition */
  auto start_decompose = std::chrono::high_resolution_clock::now();
  auto [desc_bits, asc_bits_list] = decompose(g, sub_bits, exo_bits);
  auto end_decompose = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_decompose = end_decompose - start_decompose;
  INFO("[RAI] order=" << k
        << ", time decomposing=" << elapsed_decompose.count() << "s");

  /* Stage C */
  for (auto& asc_bits : asc_bits_list) {
    rai_recursive(k + 1, asc_bits, exo_bits, g, sepset, df, test, max_cond_vars);
  }

  /* Stage D */
  std::vector<uint64_t> exo_bits_for_desc = exo_bits;
  for (auto& asc_bits : asc_bits_list) {
    for (size_t i = 0; i < exo_bits_for_desc.size(); ++i) {
      exo_bits_for_desc[i] |= asc_bits[i];
    }
  }
  rai_recursive(k + 1, desc_bits, exo_bits_for_desc, g, sepset, df, test, max_cond_vars);
}

PDAG rai(const DataframeWrapper& df,
         const CITestType& test,
         size_t max_cond_vars) {
  const size_t n = df.num_vars;
  
  PDAGwithAdjMat g(n);
  g.set_as_complete();
  Sepset sepset(n, std::vector<std::unordered_set<size_t>>(n));

  // GO!
  const size_t blocks = (n + 63) / 64;
  std::vector<uint64_t> sub_bits(blocks, ~0ULL);
  if (blocks) {
    const auto last = n % 64;
    sub_bits.back() = (last == 0) ? ~0ULL : ((1ULL << last) - 1);
  }
  std::vector<uint64_t> exo_bits(blocks, 0ULL);
  rai_recursive(0, sub_bits, exo_bits, g, sepset, df, test, max_cond_vars);

  g.orient_colliders(sepset);
  g.apply_meeks_rules();
  return g.to_pdag();
}
