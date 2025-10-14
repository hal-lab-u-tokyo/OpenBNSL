#pragma once
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "dataframe_wrapper.h"

/**
 * @ingroup base
 * @struct ContingencyTable
 * @brief Represents a contingency table for a subset of variables.
 * Key design:
 * radix_weights[i] == product of cardinalities for variables to the right of i
 * (i.e. var_ids[i+1], var_ids[i+2], ...)
 * Thus the least-significant "digit" is the last variable in var_ids.
 * make_key uses these precomputed weights directly (no extra mult loop).
 */
struct ContingencyTable {
  // column indices (ascending)
  std::vector<size_t> var_ids;

  // #distinct values per var (aligned with var_ids)
  std::vector<size_t> cardinalities;

  // suffix products per var (aligned with var_ids)
  std::vector<size_t> radix_weights;

  // linear-index => frequency
  std::unordered_map<size_t, size_t> counts;

  // for internal construction (e.g., marginalize)
  ContingencyTable() = default;

  /**
   * @brief Construct a new ContingencyTable from a subset of variables.
   * @param var_ids The column indices of the variables to include (must be
   * sorted ascending).
   * @param df The DataframeWrapper containing the data.
   */
  ContingencyTable(const std::vector<size_t>& var_ids,
                   const DataframeWrapper& df);

  /**
   * @brief Marginalize this contingency table S down to a subset T ⊆ S.
   * @param var_ids_tgt Sorted ascending subset of this->var_ids.
   * @return New ContingencyTable defined on var_ids_tgt with counts aggregated.
   * Complexity: O(nnz(S) * |T|), where nnz(S) == counts.size().
   */
  ContingencyTable marginalize_to(const std::vector<size_t>& var_ids_tgt) const;

  size_t radix_weight(size_t idx) const noexcept;

  // value of var idx encoded in key
  size_t state_of(size_t key, size_t idx) const noexcept;

  // key with var idx zeroed out
  size_t strip(size_t key, size_t idx) const noexcept;

  bool contains(size_t key) const noexcept;

  /**
   * @brief Create a linear key from a row-like object (supports operator[]).
   * Uses precomputed radix_weights for clarity and speed.
   */
  template <typename RowLike>
  size_t make_key(const RowLike& row) const noexcept {
    size_t key = 0;
    const size_t k = var_ids.size();
    for (size_t i = 0; i < k; ++i) {
      key += static_cast<size_t>(row[var_ids[i]]) * radix_weights[i];
    }
    return key;
  }
};