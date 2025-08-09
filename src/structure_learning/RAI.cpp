#include "structure_learning/RAI.h"

#include "graph/pdag_with_adjmat.h"

PDAG rai(const DataframeWrapper& df,
         const CITestType& ci_test_type,
         size_t max_cond_vars) {
  size_t n = df.num_of_vars;
  PDAGwithAdjMat g(n);

  return g.to_pdag();
}
