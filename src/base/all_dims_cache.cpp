#include "base/all_dims_cache.h"

AllDimsCache::AllDimsCache(const DataframeWrapper& df, int max_depth) : df(df) {
  if (max_depth < 0 || max_depth > (int)df.num_of_vars)
    throw std::invalid_argument("max_depth must be in [0, num_of_vars]");
  this->max_depth = max_depth;

  bit_masks.resize(df.num_of_vars);
  for (size_t i = 0; i < df.num_of_vars; ++i) {
    bit_masks[i] = (mp::cpp_int(1) << i);
  }

  root_node = std::make_unique<Node>();
  std::vector<int> freq_tbl = {(int)df.num_of_datapoints};
  std::vector<int> indices(df.num_of_datapoints);     // indices of datapoints
  std::iota(indices.begin(), indices.end(), 0);       // 0, 1, 2, ..., N - 1
  branch(root_node.get(), 0, -1, freq_tbl, indices);  // recursive
}

void AllDimsCache::branch(Node* node, const size_t curr_depth,
                          const int rightmost_var,
                          const std::vector<int>& freq_tbl,
                          const std::vector<int>& indices) {
  node->freq_tbl = freq_tbl;
  if (curr_depth >= max_depth) return;

  int num_of_branches = df.num_of_vars - rightmost_var - 1;
  node->children.resize(num_of_branches);
#pragma omp parallel for schedule(dynamic)
  for (int branch_idx = 0; branch_idx < num_of_branches; ++branch_idx) {
    const int tgt_var = rightmost_var + branch_idx + 1;
    const int num_of_values = df.val_idx2str[tgt_var].size();
    const std::vector<uint8_t>& tgt_dataref = df.data_col_idx[tgt_var];

    std::vector<int> new_freq_tbl;
    std::vector<int> new_indices(indices.size());

    int left = 0, right = 0, idx_idx = 0;
    std::vector<std::vector<int>> buckets(num_of_values);
    for (const int bucket_size : freq_tbl) {
      for (auto& bucket : buckets) bucket.clear();
      left = right;
      right += bucket_size;
      for (int j = left; j < right; ++j) {
        auto idx = indices[j];
        auto value = tgt_dataref[idx];
        buckets[value].push_back(idx);
      }

      for (const auto& bucket : buckets) {
        new_freq_tbl.push_back(bucket.size());
        for (const auto idx : bucket) {
          new_indices[idx_idx++] = idx;
        }
      }
    }

    auto new_child = std::make_unique<Node>();
    branch(new_child.get(), curr_depth + 1, tgt_var, new_freq_tbl, new_indices);
    node->children[branch_idx] = std::move(new_child);
  }
}

const std::vector<int>& AllDimsCache::get_freq_tbl(
    const std::vector<int>& variables) {
  int rightmost_var = -1;
  Node* node = root_node.get();
  for (auto tgt_var : variables) {
    auto branch_idx = tgt_var - rightmost_var - 1;
    node = node->children[branch_idx].get();
    rightmost_var = tgt_var;
  }
  return node->freq_tbl;
}