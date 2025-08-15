#pragma once
#include <variant>

/**
 * @ingroup score
 * @struct BDeu
 * @brief Structure representing the BDeu score type.
 */
struct BDeu {
  double ess;
  BDeu(double ess = 1.0) : ess(ess) {}
};

/**
 * @ingroup score
 * @struct K2
 * @brief Structure representing the K2 score type.
 */
struct K2 {};

/**
 * @ingroup score
 * @struct BIC
 * @brief Structure representing the BIC score type.
 */
struct BIC {};

/**
 * @ingroup score
 * @struct AIC
 * @brief Structure representing the AIC score type.
 */
struct AIC {};

/**
 * @ingroup score
 * @brief Type representing all possible score types.
 */
typedef std::variant<BDeu, K2, BIC, AIC> ScoreType;

/**
 * @ingroup score
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
 * @ingroup score
 * @brief Get the score type of a specific type
 * @tparam T A score type
 * @param score A score type
 * @return The score type of the specific type
 */
template <typename T>
const T& get_type(const ScoreType& score) {
  return std::get<T>(score);
}