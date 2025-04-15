#pragma once
#include <variant>

struct BDeu {
  double ess;
  BDeu(double ess = 1.0) : ess(ess) {}
};
struct K2 {};
struct BIC {};
struct AIC {};

using ScoreType = std::variant<BDeu, K2, BIC, AIC>;

/**
 * @brief Check if the score type is of a specific type
 * @tparam T A score type
 * @param score A score type
 * @return True if the score type is of the specific type; otherwise, false
 */
template <typename T>
bool is_type(const ScoreType& score) {
  return std::holds_alternative<T>(score);
}

/**
 * @brief Get the score type of a specific type
 * @tparam T A score type
 * @param score A score type
 * @return The score type of the specific type
 * @example
 * Example usage:
 * @code
 * ScoreType score = BDeu{1.0};
 * if (is_type<BDeu>(score)) {
 *     const auto& bdeu = get_type<BDeu>(score);
 *     auto ess = bdue.ess;
 *     // Process BDeu score
 * }
 * @endcode
 */
template <typename T>
const T& get_type(const ScoreType& score) {
  return std::get<T>(score);
}