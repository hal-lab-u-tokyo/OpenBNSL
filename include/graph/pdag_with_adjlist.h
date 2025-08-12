#pragma once
#include "graph/ipdag_convertible.h"
#include "graph/pdag.h"

template <bool Deterministic>
using ParentSet = std::conditional_t<Deterministic,
                                     std::set<std::size_t>,
                                     std::unordered_set<std::size_t>>;

/**
 * @ingroup graph
 * @struct PDAGwithAdjList
 * @brief PDAG implementation backed by an adjacency list.
 * @details
 * This structure uses an adjacency list to represent the PDAG, allowing for
 * efficient storage and access patterns.
 */
template <bool Deterministic>
struct PDAGwithAdjList : IPDAGConvertible {
  using ParentSetType = ParentSet<Deterministic>;

  std::size_t num_vars;
  std::vector<ParentSet<Deterministic>> parents;

  /**
   * @brief Construct a new PDAGwithAdjList with a specified number of
   * variables.
   * @param num_vars The number of variables in the PDAG.
   */
  explicit PDAGwithAdjList(std::size_t num_vars)
      : num_vars(num_vars), parents(num_vars) {}

  bool has_edge(std::size_t from, std::size_t to) const {
    return parents[from].find(to) != parents[from].end();
  }

  void add_edge(std::size_t from, std::size_t to) { parents[from].insert(to); }

  void remove_edge(std::size_t from, std::size_t to) {
    parents[from].erase(to);
  }

  bool has_path(std::size_t src, std::size_t dst) const {
    std::vector<bool> visited(this->num_vars, false);
    std::vector<std::size_t> stack{src};

    while (!stack.empty()) {
      auto v = stack.back();
      stack.pop_back();
      if (v == dst) return true;
      if (visited[v]) continue;
      visited[v] = true;

      for (auto& u : this->parents[v]) {
        if (!visited[u]) stack.push_back(u);
      }
    }
    return false;
  }

  void set_parents(std::size_t v, const ParentSetType& new_parents) {
    this->parents[v] = new_parents;
  }

  PDAG to_pdag() const {
    PDAG pdag(this->num_vars);
    for (std::size_t v = 0; v < this->num_vars; ++v) {
      for (auto p : this->parents[v]) {
        pdag.add_edge(p, v);
      }
    }
    return pdag;
  }
};