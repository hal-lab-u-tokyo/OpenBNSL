#pragma once
#include <concepts>
#include <vector>

namespace utils {

/**
 * @brief Generates the next combination of t nodes from n nodes.
 *
 * The function calculates the next combination in lexicographical order, given
 an combination represented as a bitmask.
 * The initial input should be a bitmask with rightmost `t` bits set to 1 and
 the remaining `n - t` bits set to 0.
 * For example, in case of n = 5 and t = 3, the initial combination should be
 00111 and the function will generate the following 10 combinations: 00111,
 01011, 01101, 01110, ..., 11010, 11100 (in lexicographical order).
 *
 * @param combination A reference to the current combination represented as a
 bitmask.
 * @param n The total number of nodes.
 * @return Returns `false` if all combinations have been generated; otherwise,
 `true`.

 * @example
 * Example usage:
 * @code
 * #include <boost/multiprecision/cpp_int.hpp>
 * #include "next_combmask.h"
 * namespace mp = boost::multiprecision;
 * int n = 5, t = 3;
 * mp::cpp_int combination = (mp::cpp_int(1) << t) - 1; // 00111
 * do {
 *     // Process the current combination
 * } while (next_combmask(combination, n));
 * @endcode
 *
 */
template <std::integral T>
bool next_combmask(T& comb, int n) {
  if (comb == 0) return false;
  T one = 1;
  T x = comb & -comb;                   // get the rightmost bit
  T y = comb + x;                       // move the rightmost bit to the left
  comb = (((comb & ~y) / x) >> 1) | y;  // set right bits
  return comb < (one << n);             // return false if we have done all
}

/**
 * @brief Convert a combination to a vector
 * @tparam T An integral type
 * @param combmask A combination mask represented as an integral type
 * @return A vector of indices
 */
template <std::integral T>
std::vector<size_t> combmask2vec(T combmask) {
  std::vector<size_t> vec;
  while (combmask) {
    T x = combmask & -combmask;  // get the rightmost bit
    if constexpr (sizeof(T) <= sizeof(int))
      vec.push_back(static_cast<size_t>(__builtin_ctz(x)));
    else if constexpr (sizeof(T) <= sizeof(long long))
      vec.push_back(static_cast<size_t>(__builtin_ctzll(x)));
    else
      static_assert(sizeof(T) <= sizeof(int), "T is too large");
    combmask &= combmask - 1;  // clear the rightmost bit
  }
  return vec;
}

}  // namespace utils