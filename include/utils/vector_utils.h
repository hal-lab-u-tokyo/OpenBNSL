#pragma once
#include <vector>
#include <concepts> 

namespace utils {

template <typename T>
inline std::vector<T> copy_without(const std::vector<T>& vec, const T& elem) {
  std::vector<T> res;
  res.reserve(vec.size());
  for (const auto& v : vec) {
    if (v != elem) res.push_back(v);
  }
  return res;
}

template <std::integral T>
inline std::vector<T> reindex_excluding(T x, const std::vector<T>& S) {
  std::vector<T> idx;
  idx.reserve(S.size());
  for (auto y : S) {
    idx.push_back(y < x ? y : y - 1);
  }
  return idx;
}

}  // namespace utils
