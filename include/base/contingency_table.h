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
 */
struct ContingencyTable {
  std::vector<size_t> var_ids;                // column indices (ascending)
  std::vector<size_t> cardinalities;          // #distinct values for each var
  std::vector<size_t> radix_weights;          // radix weights for each var
  std::unordered_map<size_t, size_t> counts;  // linear‑index → frequency

  /**
   * @brief Construct a new ContingencyTable from a subset of variables.
   * @param var_ids The column indices of the variables to include.
   * @param df The DataFrameWrapper containing the data.
   */
  ContingencyTable(const std::vector<size_t>& var_ids,
                   const DataframeWrapper& df);

  size_t radix_weight(size_t idx) const noexcept;
  size_t state_of(size_t key, size_t idx) const noexcept;
  size_t strip(size_t key, size_t idx) const noexcept;
  bool contains(size_t key) const noexcept;

  template <typename RowLike>
  size_t make_key(const RowLike& row) const noexcept;
};
