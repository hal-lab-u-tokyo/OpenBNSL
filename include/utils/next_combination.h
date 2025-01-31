#pragma once
#include <boost/multiprecision/cpp_int.hpp>
namespace mp = boost::multiprecision;

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
 * int n = 5, t = 3;
 * mp::cpp_int combination = (mp::cpp_int(1) << t) - 1; // 00111
 * do {
 *     // Process the current combination
 * } while (next_combination(combination, n));
 * @endcode
 *
 */
template <typename T>
bool next_combination(T& comb, int n) {
  static_assert(
      std::is_integral<T>::value || std::is_same<T, mp::cpp_int>::value,
      "T must be an integral type or cpp_int");
  if (comb == 0) return false;
  T one = 1;
  T x = comb & -comb;                   // get the rightmost bit
  T y = comb + x;                       // move the rightmost bit to the left
  comb = (((comb & ~y) / x) >> 1) | y;  // set right bits
  return comb < (one << n);             // return false if we have done all
}