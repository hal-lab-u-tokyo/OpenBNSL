#pragma once
#include <variant>

struct ChiSquare {
  double level;
  ChiSquare(double level = 0.05) : level(level) {}
};

struct GSquare {
  double level;
  GSquare(double level = 0.05) : level(level) {}
};

using CITestType = std::variant<ChiSquare, GSquare>;
/**
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
 * @brief Get the CI test type of a specific type
 * @tparam T A CI test type
 * @param ci_test A CI test type
 * @return The CI test type of the specific type
 * @example
 * Example usage:
 * @code
 * CITestType ci_test = ChiSquare{0.05};
 * if (is_type<ChiSquare>(ci_test)) {
 *     const auto& chi_square = get_type<ChiSquare>(ci_test);
 *     auto level = chi_square.level;
 *     // Process ChiSquare CI test
 * }
 * @endcode
 */
template <typename T>
const T& get_type(const CITestType& ci_test) {
  return std::get<T>(ci_test);
}