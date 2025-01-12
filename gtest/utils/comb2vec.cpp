#include "utils/comb2vec.h"

#include <gtest/gtest.h>

TEST(Comb2VecTest, correctness) {
  // Test case 1: Empty combination
  {
    unsigned int comb = 0;  // No bits set
    std::vector<int> result = comb2vec(comb);
    ASSERT_TRUE(result.empty());
  }

  // Test case 2: Single bit
  {
    unsigned int comb = 1;  // 00001
    std::vector<int> result = comb2vec(comb);
    ASSERT_EQ(result, std::vector<int>({0}));
  }

  // Test case 3: Multiple bits
  {
    unsigned int comb = (1U << 3) | (1U << 1);  // 01010
    std::vector<int> result = comb2vec(comb);
    ASSERT_EQ(result, std::vector<int>({1, 3}));
  }

  // Test case 4: All bits set in a small range
  {
    unsigned int comb = (1U << 5) - 1;  // 11111 (bits 0-4)
    std::vector<int> result = comb2vec(comb);
    ASSERT_EQ(result, std::vector<int>({0, 1, 2, 3, 4}));
  }

  // Test case 5: Sparse bits with 64-bit input
  {
    unsigned long long comb =
        (1ULL << 50) | (1ULL << 30);  // Large numbers with sparse bits
    std::vector<int> result = comb2vec(comb);
    ASSERT_EQ(result, std::vector<int>({30, 50}));
  }

  // Test case 6: All bits set for a 64-bit integer
  {
    unsigned long long comb = ~0ULL;  // All bits set
    std::vector<int> result = comb2vec(comb);
    std::vector<int> expected(64);
    for (int i = 0; i < 64; ++i) expected[i] = i;
    ASSERT_EQ(result, expected);
  }

  // Test case 7: Unsupported type (should trigger static_assert)
  // This cannot be directly tested at runtime, but static_assert will ensure
  // compile-time correctness.
}
