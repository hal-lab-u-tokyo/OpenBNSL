#include "util.h"

#include <gmpxx.h>
#include <gtest/gtest.h>

int comb_cnt(int n, int m) {
  mpz_class comb = (mpz_class(1) << m) - 1;
  int cnt = 0;
  do {
    cnt++;
  } while (next_combination(comb, n));
  return cnt;
}

TEST(UtilTest, correctness) {
  ASSERT_EQ(comb_cnt(5, 0), 1);
  ASSERT_EQ(comb_cnt(5, 3), 10);
  ASSERT_EQ(comb_cnt(10, 5), 252);
  ASSERT_EQ(comb_cnt(20, 10), 184756);
}