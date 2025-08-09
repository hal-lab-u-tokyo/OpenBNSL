#include <gtest/gtest.h>

#include <boost/multiprecision/cpp_int.hpp>

#include "utils/next_combmask.h"
namespace mp = boost::multiprecision;

template <typename T>
int comb_cnt(int n, int m) {
  T comb = (T(1) << m) - 1;  // 初期ビットマスク（右端 m ビットが 1）
  int cnt = 0;
  do {
    cnt++;
  } while (next_combmask(comb, n));
  return cnt;
}

// Test case 1: 32 bit integer
TEST(NextCombinationTest, StandardInt) {
  ASSERT_EQ(comb_cnt<unsigned int>(5, 3), 10);
  ASSERT_EQ(comb_cnt<unsigned int>(10, 3), 120);
  ASSERT_EQ(comb_cnt<unsigned int>(31, 3), 4495);
}

// Test case 2: 64 bit integer
TEST(NextCombinationTest, LongLongInt) {
  ASSERT_EQ(comb_cnt<unsigned long long>(5, 3), 10);
  ASSERT_EQ(comb_cnt<unsigned long long>(10, 3), 120);
  ASSERT_EQ(comb_cnt<unsigned long long>(31, 3), 4495);
  ASSERT_EQ(comb_cnt<unsigned long long>(63, 3), 39711);
}

// Test case 3: Multiprecision integer
TEST(NextCombinationTest, MultiprecisionCppInt) {
  ASSERT_EQ(comb_cnt<mp::cpp_int>(5, 3), 10);
  ASSERT_EQ(comb_cnt<mp::cpp_int>(10, 3), 120);
  ASSERT_EQ(comb_cnt<mp::cpp_int>(31, 3), 4495);
  ASSERT_EQ(comb_cnt<mp::cpp_int>(63, 3), 39711);
  ASSERT_EQ(comb_cnt<mp::cpp_int>(127, 3), 333375);
}

// Test case 4: Boundary conditions
TEST(NextCombinationTest, BoundaryConditions) {
  // n = 0, m = 0
  ASSERT_EQ(comb_cnt<unsigned int>(0, 0), 1);
  ASSERT_EQ(comb_cnt<unsigned long long>(0, 0), 1);
  ASSERT_EQ(comb_cnt<mp::cpp_int>(0, 0), 1);

  // n > 0, m = 0
  ASSERT_EQ(comb_cnt<unsigned int>(31, 0), 1);
  ASSERT_EQ(comb_cnt<unsigned long long>(63, 0), 1);
  ASSERT_EQ(comb_cnt<mp::cpp_int>(127, 0), 1);

  // n > 0, m = n
  ASSERT_EQ(comb_cnt<unsigned int>(31, 31), 1);
  ASSERT_EQ(comb_cnt<unsigned long long>(63, 63), 1);
  ASSERT_EQ(comb_cnt<mp::cpp_int>(127, 127), 1);
}