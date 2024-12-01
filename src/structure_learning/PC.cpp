#include "structure_learning/PC.h"

#include <gmpxx.h>

#include <vector>

#include "base/PDAG.h"
#include "util.h"

// input: data: np.ndarray,  shape: (n: number of variables, d: number of
// samples) output: leard PDAG
PDAGwithAdjMat PC(const py::array_t<int> &data, int max_order) {
  auto buf = data.request();
  int *ptr = (int *)buf.ptr;
  size_t n = buf.shape[0];
  size_t d = buf.shape[1];

  std::vector<int> variables;
  for (int i = 0; i < n; i++) {
    variables.push_back(i);
  }
  PDAGwithAdjMat pdag(n);

  // Phase 1: Skeleton
  int order = 0;
  while (order < n - 2) {
    for (int x_i = 0; x_i < n; x_i++) {
      for (int x_j = x_i + 1; x_j < n; x_j++) {
        std::vector<int> common_neighbors;
        for (auto x_k : variables) {
          if (x_k == x_i || x_k == x_j) {
            continue;
          }
          if (pdag.has_edge(x_i, x_k) && pdag.has_edge(x_j, x_k)) {
            common_neighbors.push_back(x_k);
          }
        }

        int cn_size = common_neighbors.size();
        mpz_class comb_set = (mpz_class(1) << order) - 1;
        do {
          // if ( CItest(x_i, x_j, comb_set, data) == 1 ) pdag.remove_edge(x_i,
          // x_j);
        } while (next_combination(comb_set, cn_size));
      }
    }
    if (order == max_order) break;
    order++;
  }

  // Phase 2: Orientation

  // detect all v-structures
  // for G でX <-> Y <-> Z が形成するすべての3ノードの組み合わせに対して
  // もし，ＸとＺが隣接せず，CItest(X, Z, Y, data) == 0 (従属) ならば，
  // X -> Y <- Z が確定する
  // 無向辺は，X <-> Y <-> Z で表現しているので，X <- Y と Y -> をremove_edge

  // the following should be implemented in PDAG class

  // while (方向付けられる辺がある){
  //     if X -> Y <-> Z かつ，XとZが隣接していないならば，
  //     V構造として検出されていないため，X -> Y -> Z が確定，Y <- Z
  //     をremove_edge

  //     if X -> Y -> Z かつ，X <-> Z ならば，
  //     巡回rを作らないように，X -> Z が確定，X <- Z をremove_edge

  //     if X -> Y <- Z かつ，XとZが隣接していないならば，

  // }

  return pdag;
}
