#pragma once
#include <cstddef>
#include <vector>

#include "dataframe_wrapper.h"

struct ContingencyTable {
  std::vector<size_t> var_ids;
  std::vector<size_t> cardinalities;
  std::vector<size_t> counts;
};

ContingencyTable buildContingencyTable(const std::vector<size_t>& var_ids,
                                       const DataframeWrapper& df);
