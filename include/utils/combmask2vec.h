#pragma once
#include <cstddef>
#include <vector>

/**
 * @brief Convert a combination to a vector
 * @tparam T An integral type
 * @param combmask A combination mask represented as an integral type
 * @return A vector of indices
 */
template <typename T>
std::vector<size_t> combmask2vec(T combmask) {
  static_assert(std::is_integral<T>::value, "T must be an integral type");
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
