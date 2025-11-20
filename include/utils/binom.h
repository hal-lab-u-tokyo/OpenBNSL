#pragma once
#include <vector>

namespace utils {

/**
 * @brief Build a binomial coefficient table up to max_n and max_k.
 * @param max_n Maximum n value.
 * @param max_k Maximum k value.
 * @return A 2D vector where the entry at [n][k] is C(n, k).
 */
inline std::vector<std::vector<size_t>> build_binom_table(size_t max_n,
                                                          size_t max_k) {
  std::vector<std::vector<size_t>> res(max_n + 1);
  for (size_t n = 0; n <= max_n; ++n) {
    res[n].resize(max_k + 1);
    res[n][0] = 1;
    for (size_t k = 1; k <= std::min(n, max_k); ++k) {
      res[n][k] = res[n - 1][k - 1] + res[n - 1][k];
    }
  }
  return res;
}

}  // namespace utils