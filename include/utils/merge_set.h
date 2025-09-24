#pragma once
#include <algorithm>
#include <vector>

template <typename T>
static std::vector<T> merge_set(const std::vector<T>& a,
                                const std::vector<T>& b) {
  std::vector<T> r;
  r.reserve(a.size() + b.size());
  r.insert(r.end(), a.begin(), a.end());
  r.insert(r.end(), b.begin(), b.end());
  std::sort(r.begin(), r.end());
  r.erase(std::unique(r.begin(), r.end()), r.end());
  return r;
}