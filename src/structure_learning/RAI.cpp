using namespace std;
//#include "structure_learning/RAI.h"

//#include <gmpxx.h>
#include <string>
#include <vector>

//#include "base/PDAG.h"
//#include "util.h"

// input: data: np.ndarray,  shape: (n: number of variables, d: number of
// samples) output: leard PDAG
vector<vector<bool>> RAI(const vector<vector<string>> &data, vector<vector<bool>> &Gall, vector<int> Gs, vector<int> Gex, int N, string ci_test, float ESS, bool first) {
  int n_node = data.at(0).size();
  if (first) {// initialize Gall, Gs, Gex
      vector<vector<bool>> gall(n_node, vector<bool>(n_node, 1));
      Gall = gall; //これで参照元を変更できる？ Gall is complete graph
      vector<int> gs(n_node, 1);
      for (int i = 0; i < n_node; i++) {
          gs.at(i) = i; // Gs contains all nodes 0 ~ n_node-1
      }
      vector<int> gex;
      Gex = gex; // Gex is empty
      first = false;
  }


  return Gall;
}
