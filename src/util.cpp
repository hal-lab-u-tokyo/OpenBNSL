#include <gmpxx.h>

// get all combinations of t nodes from n nodes
// initial input: 000...011...1 (n-t 0s and t 1s)
bool next_combination(mpz_class& combination, int n) {
    if (combination == 0) return false; 
    mpz_class x = combination & -combination;           // get the rightmost bit
    mpz_class y = combination + x;                      // move the rightmost bit to the left
    combination = (((combination & ~y) / x) >> 1) | y;  // set right bits
    return combination < mpz_class(1) << n;             // return false if we have done all combinations
}

// 00111 => {2, 1, 0}
// 01011 => {3, 1, 0}
// 01101 => {3, 2, 0}
// 01110 => {3, 2, 1}
// ...
// 11010 => {4, 3, 1}
// 11100 => {4, 3, 2}