#pragma once
#include <cstddef>
#include <set>
#include <vector>

/**
 * @ingroup graph
 * @struct PDAG
 * @brief Interface for converting internal graph structures to a unified PDAG.
 * @details
 * Structure learning algorithms may use different graph representations.
 * This interface ensures that all such graphs can be converted into a common
 * sparse PDAG form for downstream use (e.g., evaluation, comparison, or Python
 * interop).
 */
struct PDAG {
  /* Data members */
  std::size_t num_vars;
  std::vector<std::set<size_t>> parents;  // parents[v] = {u | u -> v}

  /* Lifecycle */
  /**
   * @brief Construct a new PDAG with a specified number of variables.
   * @param num_vars The number of variables in the PDAG.
   */
  PDAG(std::size_t num_vars) : num_vars(num_vars), parents(num_vars) {}

  bool has_edge(std::size_t from, std::size_t to) const {
    return parents[to].find(from) != parents[to].end();
  }

  void add_edge(std::size_t from, std::size_t to) { parents[to].insert(from); }

  void remove_edge(std::size_t from, std::size_t to) {
    parents[to].erase(from);
  }
};
