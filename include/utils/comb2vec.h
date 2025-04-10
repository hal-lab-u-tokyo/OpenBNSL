#pragma once
#include <cstddef>
#include <vector>

/**
 * @brief Convert a combination to a vector
 * @tparam T An integral type
 * @param comb A combination
 * @return A vector of indices
 */
template <typename T>
std::vector<size_t> comb2vec(T comb) {
  static_assert(std::is_integral<T>::value, "T must be an integral type");
  std::vector<size_t> vec;
  while (comb) {
    T x = comb & -comb;  // get the rightmost bit
    if constexpr (sizeof(T) <= sizeof(int))
      vec.push_back(static_cast<size_t>(__builtin_ctz(x)));
    else if constexpr (sizeof(T) <= sizeof(long long))
      vec.push_back(static_cast<size_t>(__builtin_ctzll(x)));
    else
      static_assert(sizeof(T) <= sizeof(int), "T is too large");
    comb &= comb - 1;  // clear the rightmost bit
  }
  return vec;
}
