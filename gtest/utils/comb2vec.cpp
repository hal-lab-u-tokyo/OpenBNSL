#include "utils/comb2vec.h"

#include <gtest/gtest.h>

#include <vector>

TEST(Comb2VecTest, correctness) {
  // Test case 1: Empty combination
  {
    unsigned int comb = 0;  // No bits set
    std::vector<size_t> result = comb2vec(comb);
    ASSERT_TRUE(result.empty());
  }

  // Test case 2: Single bit
  {
    unsigned int comb = 1;  // 00001
    std::vector<size_t> result = comb2vec(comb);
    std::vector<size_t> expected = {0};
    ASSERT_EQ(result, expected);
  }

  // Test case 3: Multiple bits
  {
    unsigned int comb = (1U << 3) | (1U << 1);  // 0b01010, bits 1 and 3 set
    std::vector<size_t> result = comb2vec(comb);
    std::vector<size_t> expected = {1, 3};
    ASSERT_EQ(result, expected);
  }

  // Test case 4: All bits set in a small range
  {
    unsigned int comb = (1U << 5) - 1;  // 0b11111 (bits 0-4)
    std::vector<size_t> result = comb2vec(comb);
    std::vector<size_t> expected = {0, 1, 2, 3, 4};
    ASSERT_EQ(result, expected);
  }

  // Test case 5: Sparse bits with 64-bit input
  {
    unsigned long long comb =
        (1ULL << 50) | (1ULL << 30);  // Sparse bits: positions 30 and 50
    std::vector<size_t> result = comb2vec(comb);
    std::vector<size_t> expected = {30, 50};
    ASSERT_EQ(result, expected);
  }

  // Test case 6: All bits set for a 64-bit integer
  {
    unsigned long long comb = ~0ULL;  // All 64 bits set
    std::vector<size_t> result = comb2vec(comb);
    std::vector<size_t> expected;
    for (size_t i = 0; i < 64; ++i) expected.push_back(i);
    ASSERT_EQ(result, expected);
  }
}
