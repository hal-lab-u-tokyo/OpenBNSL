#pragma once
#include <cstddef>
#include <vector>

#include "dataframe_wrapper.h"

struct ContingencyTable {
  std::vector<size_t> vars;
  std::vector<size_t> cardinalities;
  std::vector<size_t> counts;
};

ContingencyTable buildContingencyTable(const std::vector<size_t>& vars,
                                       const DataframeWrapper& df);
