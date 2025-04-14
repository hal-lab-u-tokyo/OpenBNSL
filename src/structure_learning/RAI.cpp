#include "structure_learning/RAI.h"

PDAG RAI(const DataframeWrapper& df, const CITestType& ci_test_type) {
  size_t n = df.num_of_vars;
  PDAG g(n);
  g.complete_graph();

  return g;
}
