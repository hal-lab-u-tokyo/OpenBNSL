#pragma once
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace utils {

template <typename T>
void enum_combvec_rec(const std::vector<T>& items,
                      size_t k,
                      size_t idx,
                      std::vector<T>& cur,
                      std::vector<std::vector<T>>& out) {
  if (cur.size() == k) {
    out.push_back(cur);
    return;
  }
  for (size_t i = idx; i < items.size(); ++i) {
    cur.push_back(items[i]);
    enum_combvec_rec(items, k, i + 1, cur, out);
    cur.pop_back();
  }
}

/**
 * @brief Generate all combinations of size k from a vector of items
 * @tparam T The type of items in the vector
 * @param items A vector of items to generate combinations from
 * @return A vector of vectors, where each inner vector is a combination of size
 *
 * @example
 * Example usage:
 * @code
 * std::vector<size_t> items = {0, 1, 2, 3};
 * auto combs = utils::gen_combvecs(items, 2);
 * // => {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}}
 * @endcode
 */
template <typename T>
std::vector<std::vector<T>> gen_combvecs(const std::vector<T>& items,
                                         size_t k) {
  if (!std::is_sorted(items.begin(), items.end()))
    throw std::invalid_argument("items must be sorted in ascending order");
  std::vector<std::vector<T>> res;
  if (k == 0) {
    res.push_back({});
    return res;
  }
  std::vector<T> cur;
  enum_combvec_rec(items, k, 0, cur, res);
  return res;
}

/**
 * @brief Compute the lexicographic rank of a combination vector.
 * @tparam T Integral index type (e.g., size_t, int, uint32_t)
 * @param S Combination vector (must be sorted)
 * @param nCk_tbl Precomputed nCk table (as from nCk_tbl)
 * @return Lexicographic rank of S
 */
template <std::integral T>
size_t combvec_rank(const std::vector<T>& S,
                    const std::vector<std::vector<size_t>>& nCk_tbl) {
  size_t rank = 0;
  for (size_t i = 0; i < S.size(); ++i) {
    rank += nCk_tbl[S[i]][i + 1];
  }
  return rank;
}

}  // namespace utils