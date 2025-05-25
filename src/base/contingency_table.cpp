#include "base/contingency_table.h"

ContingencyTable buildContingencyTable(const std::vector<size_t>& var_ids,
                                       const DataframeWrapper& df) {
  if (!std::is_sorted(var_ids.begin(), var_ids.end())) {
    throw std::invalid_argument("Input variable indices must be sorted.");
  }

  ContingencyTable ct;
  ct.var_ids = var_ids;

  size_t totalCells = 1;

  for (size_t v : var_ids) {
    size_t nvals = df.num_of_values[v];
    ct.cardinalities.push_back(nvals);
    totalCells *= nvals;
  }
  ct.counts.resize(totalCells, 0);

  for (size_t i = 0; i < df.num_of_datapoints; i++) {
    size_t idx = 0;
    size_t multiplier = 1;
    for (int j = static_cast<int>(var_ids.size()) - 1; j >= 0; j--) {
      size_t v = var_ids[j];
      idx += df.data_row_major[i][v] * multiplier;
      multiplier *= ct.cardinalities[j];
    }
    ct.counts[idx]++;
  }
  return ct;
}
