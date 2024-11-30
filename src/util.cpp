#include <gmpxx.h>

bool next_combination(mpz_class& combination, int n) {
    if (combination == 0) return false; 
    mpz_class x = combination & -combination;           // get the rightmost bit
    mpz_class y = combination + x;                      // move the rightmost bit to the left
    combination = (((combination & ~y) / x) >> 1) | y;  // set right bits
    return combination < mpz_class(1) << n;             // return false if we have done all combinations
}