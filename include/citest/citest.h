#pragma once

#include <cstddef>
#include <vector>

#include "base/contingency_table.h"
#include "citest/citest_type.h"

/**
 * @ingroup citest
 * @brief Conditional–independence test for discrete variables.
 * @details
 * The contingency table @p ct is stored sparsely, i.e. cells with observed
 * count 0 do not appear.  pgmpy/Scipy, however, still accounts for such cells
 * if the expected frequency is non‑zero.  We therefore generate those cells
 * on‑the‑fly: for each Z‑slice we iterate over the Cartesian product of the
 * actually observed X–states and Y–states and fill in obs = 0 where necessary.
 * @param x The first variable
 * @param y The second variable
 * @param sepset_candidate The candidate separator set
 * @param ct The contingency table
 * @param ci_test_type The CI test type
 * @return True if the conditional independence test passes (p_value >= alpha),
 *         false otherwise.
 */
bool citest(std::size_t x,
            std::size_t y,
            const std::vector<std::size_t>& sepset_candidate,
            const ContingencyTable& ct,
            const CITestType& ci_test_type);