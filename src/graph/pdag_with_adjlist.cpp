#include "graph/pdag_with_adjlist.h"
#include <queue>

template <bool Deterministic>
PDAGwithAdjList<Deterministic>::PDAGwithAdjList(std::size_t n)
    : num_vars(n), parents(n) {}

template <bool Deterministic>
bool PDAGwithAdjList<Deterministic>::has_edge(std::size_t from, std::size_t to) const {
  const auto& ps = parents[to];
  return ps.find(from) != ps.end();
}

template <bool Deterministic>
void PDAGwithAdjList<Deterministic>::add_edge(std::size_t from, std::size_t to) {
  parents[to].insert(from);
}

template <bool Deterministic>
void PDAGwithAdjList<Deterministic>::remove_edge(std::size_t from, std::size_t to) {
  parents[to].erase(from);
}

template <bool Deterministic>
bool PDAGwithAdjList<Deterministic>::has_path(std::size_t src, std::size_t dst) const {
  if (src == dst) return true;

  std::vector<bool> visited(num_vars, false);
  std::vector<std::size_t> stack;
  stack.push_back(dst);

  while (!stack.empty()) {
    auto v = stack.back();
    stack.pop_back();
    if (visited[v]) continue;
    visited[v] = true;

    for (auto p : parents[v]) {
      if (p == src) return true;
      if (!visited[p]) stack.push_back(p);
    }
  }
  return false;
}

template <bool Deterministic>
void PDAGwithAdjList<Deterministic>::set_parents(std::size_t v, const ParentSetType& new_parents) {
  parents[v] = new_parents;
}

template <bool Deterministic>
PDAG PDAGwithAdjList<Deterministic>::to_pdag() const {
  PDAG pdag(num_vars);
  for (std::size_t v = 0; v < num_vars; ++v) {
    for (auto p : parents[v]) {
      pdag.add_edge(p, v);
    }
  }
  return pdag;
}

template struct PDAGwithAdjList<true>;
template struct PDAGwithAdjList<false>;