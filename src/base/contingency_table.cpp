#include "base/contingency_table.h"

#include <algorithm>
#include <stdexcept>

ContingencyTable::ContingencyTable(const std::vector<size_t>& var_ids,
                                   const DataframeWrapper& df)
    : var_ids(var_ids) {
  if (!std::is_sorted(var_ids.begin(), var_ids.end()))
    throw std::invalid_argument("var_ids must be sorted");

  cardinalities.reserve(var_ids.size());
  radix_weights.resize(var_ids.size());

  size_t mult = 1;
  for (int i = static_cast<int>(var_ids.size()) - 1; i >= 0; --i) {
    size_t v = var_ids[i];
    cardinalities.push_back(df.num_of_values[v]);
    radix_weights[i] = mult;
    mult *= df.num_of_values[v];
  }
  std::reverse(cardinalities.begin(), cardinalities.end());

  for (size_t i = 0; i < df.num_of_datapoints; ++i) {
    size_t key = make_key(df.data_row_major[i]);
    ++counts[key];
  }
}

/* ---------- key helpers --------------------------------------------- */

size_t ContingencyTable::radix_weight(size_t idx) const noexcept {
  return radix_weights[idx];
}

// value of variable `idx` encoded in `key`
size_t ContingencyTable::state_of(size_t key, size_t idx) const noexcept {
  size_t w = radix_weight(idx);
  return (key / w) % cardinalities[idx];
}

// key after zeroing‐out variable `idx`
size_t ContingencyTable::strip(size_t key, size_t idx) const noexcept {
  return key - state_of(key, idx) * radix_weight(idx);
}

bool ContingencyTable::contains(size_t key) const noexcept {
  return counts.find(key) != counts.end();
}

template <typename RowLike>
size_t ContingencyTable::make_key(const RowLike& row) const noexcept {
  size_t key = 0, mult = 1;
  for (int i = static_cast<int>(var_ids.size()) - 1; i >= 0; --i) {
    size_t v = var_ids[i];
    key += static_cast<size_t>(row[v]) * mult;
    mult *= cardinalities[i];
  }
  return key;
}