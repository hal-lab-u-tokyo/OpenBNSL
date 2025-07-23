#pragma once
#include <cstddef>
#include <set>
#include <vector>

/**
 * @brief Common interface for returning Partially Directed Acyclic Graphs
 * (PDAGs) from structure learning algorithms. Since structure learning
 * algorithms typically output sparse graphs, adjacency lists are used for
 * representation.
 */
struct PDAG {
  std::size_t num_vars;
  std::vector<std::set<size_t>> parents;  // parents[v] = {u | u -> v}

  explicit PDAG(std::size_t num_vars) : num_vars(num_vars), parents(num_vars) {}

  bool has_edge(std::size_t from, std::size_t to) const {
    return parents[from].find(to) != parents[from].end();
  }

  void add_edge(std::size_t from, std::size_t to) { parents[from].insert(to); }

  void remove_edge(std::size_t from, std::size_t to) {
    parents[from].erase(to);
  }
};
