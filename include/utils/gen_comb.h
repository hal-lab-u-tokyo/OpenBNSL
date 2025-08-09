#pragma once
#include <vector>

template <typename T>
void enum_comb_rec(const std::vector<T>& items,
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
    enum_comb_rec(items, k, i + 1, cur, out);
    cur.pop_back();
  }
}

/**
 * @brief Generate all combinations of size k from a vector of items
 * @tparam T The type of items in the vector
 * @param items A vector of items to generate combinations from
 * @return A vector of vectors, where each inner vector is a combination of size
 * k
 */
template <typename T>
std::vector<std::vector<T>> gen_combs(const std::vector<T>& items, size_t k) {
  std::vector<std::vector<T>> res;
  if (k == 0) {
    res.push_back({});
    return res;
  }
  std::vector<T> cur;
  enum_comb_rec(items, k, 0, cur, res);
  return res;
}