#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <vector>

#include "base/PDAG.h"
namespace py = pybind11;

// I. Tsamardinos, L.E. Brown, and C.F. Aliferis, “The max-min hill-climbing
// Bayesian network structure learning algorithm,” Mach. Learn., vol.65, no.1,
// pp.31–78, 2006. 9.1.4.  Measures of
// performanceにおいて，構造ハミング距離が評価基準として提案されている
// ・2つのPDAG間の構造ハミング距離は，無向辺の追加または削除、辺の向きの追加、削除、または反転の操作を必要とする回数で，アルゴリズム4に詳細がある．
// ・DAGを返すアルゴリズムの場合、この指標を計算する前にChickering,
// 1995のアルゴリズムを使用して対応するPDAGに変換する．
// ・SHDをDAGではなくPDAGで定義する理由は、統計的に区別できない構造の違いに対してペナルティを課さないようにするため
// SHD >= ME + EE, ME: missing edges, EE: extra edges
