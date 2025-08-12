#pragma once
#include <variant>

#include "graph/pdag.h"

/**
 * @ingroup citest
 * @struct ChiSquare
 * @brief CItest type for Chi-Square tests.
 */
struct ChiSquare {
  double level;
  ChiSquare(double level = 0.05) : level(level) {}
};

/**
 * @ingroup citest
 * @struct GSquare
 * @brief CItest type for G-Square tests.
 */
struct GSquare {
  double level;
  GSquare(double level = 0.05) : level(level) {}
};

/**
 * @ingroup citest
 * @struct OracleGraph
 * @brief CItest type for Oracle Graph.
 */
struct OracleGraph {
  PDAG graph;
  explicit OracleGraph(const PDAG& g) : graph(g) {}
};

/**
 * @ingroup citest
 * @typedef CITestType
 * @brief CItest type for various CI tests.
 */
typedef std::variant<ChiSquare, GSquare, OracleGraph> CITestType;

/**
 * @ingroup citest
 * @brief Check if the CI test type is of a specific type
 * @tparam T A CI test type
 * @param ci_test A CI test type
 * @return True if the CI test type is of the specific type; otherwise, false
 */
template <typename T>
bool is_type(const CITestType& ci_test) {
  return std::holds_alternative<T>(ci_test);
}

/**
 * @ingroup citest
 * @brief Get the CI test type of a specific type
 * @tparam T A CI test type
 * @param ci_test A CI test type
 * @return The CI test type of the specific type
 */
template <typename T>
const T& get_type(const CITestType& ci_test) {
  return std::get<T>(ci_test);
}