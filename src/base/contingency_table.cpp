#include "base/contingency_table.h"

#include <algorithm>
#include <stdexcept>

ContingencyTable::ContingencyTable(const std::vector<size_t>& var_ids,
                                   const DataframeWrapper& df)
    : var_ids(var_ids) {
  if (!std::is_sorted(this->var_ids.begin(), this->var_ids.end()))
    throw std::invalid_argument("var_ids must be sorted");

  const size_t k = this->var_ids.size();
  cardinalities.resize(k);
  radix_weights.resize(k);

  // 1) Fill cardinalities in-place (no reverse needed).
  for (size_t i = 0; i < k; ++i) {
    const size_t v = this->var_ids[i];
    cardinalities[i] = df.num_of_values[v];
  }

  // 2) Build suffix prods for radix_weights: weight of the last var is 1.
  size_t mult = 1;
  for (size_t i = k; i-- > 0;) {
    radix_weights[i] = mult;
    mult *= cardinalities[i];
  }

  // 3) Aggregate counts by linearized key.
  for (size_t i = 0; i < df.num_datapoints; ++i) {
    const size_t key = make_key(df.data_row_major[i]);
    ++counts[key];
  }
}

ContingencyTable ContingencyTable::marginalize_to(
    const std::vector<size_t>& var_ids_tgt) const {
  if (!std::is_sorted(var_ids_tgt.begin(), var_ids_tgt.end()))
    throw std::invalid_argument("var_ids_tgt must be sorted");

  const size_t num_vars_src = var_ids.size();
  const size_t num_vars_tgt = var_ids_tgt.size();
  std::vector<size_t> pos_in_src;  // idx_tgt -> idx_src
  pos_in_src.reserve(num_vars_tgt);
  size_t idx_tgt = 0;  // for var_ids_tgt
  size_t idx_src = 0;  // for var_ids
  while (idx_tgt < num_vars_tgt && idx_src < num_vars_src) {
    if (var_ids_tgt[idx_tgt] == var_ids[idx_src]) {
      pos_in_src.push_back(idx_src);
      ++idx_tgt;
      ++idx_src;
    } else if (var_ids_tgt[idx_tgt] > var_ids[idx_src]) {
      ++idx_src;
    } else {
      throw std::invalid_argument("var_ids_tgt must be a subset of var_ids");
    }
  }
  if (idx_tgt != num_vars_tgt) {
    throw std::invalid_argument("var_ids_tgt must be a subset of var_ids");
  }

  // Build ct_tgt
  ContingencyTable ct_tgt;
  ct_tgt.var_ids = var_ids_tgt;

  // cardinalities for ct_tgt
  ct_tgt.cardinalities.resize(num_vars_tgt);
  for (size_t t = 0; t < num_vars_tgt; ++t) {
    size_t card = cardinalities[pos_in_src[t]];
    if (card == 0) throw std::runtime_error("zero-cardinality variable");
    ct_tgt.cardinalities[t] = card;
  }

  // radix_weights for ct_tgt
  ct_tgt.radix_weights.resize(num_vars_tgt);
  size_t mult = 1;
  for (size_t t = num_vars_tgt; t-- > 0;) {
    ct_tgt.radix_weights[t] = mult;
    mult *= ct_tgt.cardinalities[t];
  }

  // Reserve counts to reduce rehash: min(nnz(S), Π card(T))
  size_t num_non_zero = counts.size();  // nnz(S)
  size_t prod = 1;
  for (size_t card : ct_tgt.cardinalities) {
    // prod * card > nnz(S) <=> prod > nnz(S) / card
    if (prod > num_non_zero / card) break;
    prod *= card;
  }
  const size_t cap = std::min(num_non_zero, prod);
  ct_tgt.counts.reserve(cap);

  // Preload arrays for tighter inner loop
  std::vector<size_t> rdxw_src(num_vars_tgt);
  std::vector<size_t> card_src(num_vars_tgt);
  std::vector<size_t> rdxw_tgt(num_vars_tgt);
  for (size_t idx_tgt = 0; idx_tgt < num_vars_tgt; ++idx_tgt) {
    const size_t pos = pos_in_src[idx_tgt];
    rdxw_src[idx_tgt] = radix_weights[pos];
    card_src[idx_tgt] = cardinalities[pos];
    rdxw_tgt[idx_tgt] = ct_tgt.radix_weights[idx_tgt];
  }

  // Fold S into T in one pass over nnz(S).
  for (const auto& [key_src, cnt_src] : counts) {
    size_t key_tgt = 0;
    for (size_t idx_tgt = 0; idx_tgt < num_vars_tgt; ++idx_tgt) {
      const size_t state = (key_src / rdxw_src[idx_tgt]) % card_src[idx_tgt];
      key_tgt += state * rdxw_tgt[idx_tgt];
    }
    ct_tgt.counts[key_tgt] += cnt_src;
  }

  return ct_tgt;
}

/* ---------- key helpers --------------------------------------------- */

size_t ContingencyTable::radix_weight(size_t idx) const noexcept {
  return radix_weights[idx];
}

// value of variable idx encoded in key
size_t ContingencyTable::state_of(size_t key, size_t idx) const noexcept {
  const size_t w = radix_weight(idx);
  return (key / w) % cardinalities[idx];
}

// key after zeroing-out variable idx
size_t ContingencyTable::strip(size_t key, size_t idx) const noexcept {
  return key - state_of(key, idx) * radix_weight(idx);
}

bool ContingencyTable::contains(size_t key) const noexcept {
  return counts.find(key) != counts.end();
}
