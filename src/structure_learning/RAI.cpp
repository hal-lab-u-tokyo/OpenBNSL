using namespace std;

#include <string>
#include <set>
#include <vector>
//the following includes are for permutation and combination algorithms
#include <algorithm>
#include <functional>
#include <iostream>


//#include "structure_learning/RAI.h"
//#include <gmpxx.h>
//#include "base/PDAG.h"
//#include "util.h"

// input: data: np.ndarray,  shape: (n: number of variables, d: number of
// samples) output: leard PDAG
// Gall.at(i).at(j)==1 means there is an edge i -> j
// c++ 17 を仮定]

void recursive_comb(int *indexes, int s, int rest, function<void(int *)> f) {
  if (rest == 0) {
    f(indexes);
  } else {
    if (s < 0) return;
    recursive_comb(indexes, s - 1, rest, f);
    indexes[rest - 1] = s;
    recursive_comb(indexes, s - 1, rest - 1, f);
  }
}

void foreach_comb(int n, int k, function<void(int *)> f) {
  int indexes[k];
  recursive_comb(indexes, n - 1, k, f);
} // foreach_comb(n, k, [](int *indexes) {indexes holds the combination nCk})


struct PDAG {
  vector<vector<bool>> g;

  vector<int> successors(int i) {
    // return the list of successors of node i
    vector<int> succ;
    for (int j = 0; j < g.size(); j++) {
      if (g.at(i).at(j)) {
        succ.push_back(j);
      }
    }
    return succ;
  }

  vector<int> predecessors(int i) {
    // return the list of predecessors of node i
    vector<int> pred;
    for (int j = 0; j < g.size(); j++) {
      if (g.at(j).at(i)) {
        pred.push_back(j);
      }
    }
    return pred;
  }

  vector<int> undirected_neighbors(int i) {
    // return the list of undirected neighbors of node i
    vector<int> neigh;
    for (int j = 0; j < g.size(); j++) {
      if (g.at(i).at(j) && g.at(j).at(i)) {
        neigh.push_back(j);
      }
    }
    return neigh;
  }

  void remove_oriented_edge(int i, int j) {
    // remove the edge i -> j
    g.at(i).at(j) = false;
  }

  void remove_edge(int i, int j) {
    // remove edge between i and j
    g.at(i).at(j) = false;
    g.at(j).at(i) = false;
  }

  void add_edge(int i, int j) {
    g.at(i).at(j) = true;
  }

  bool has_edge(int i, int j) {
    return g.at(i).at(j);
  }

  bool has_directed_edge(int i, int j) {
    if (g.at(i).at(j) && !g.at(j).at(i)) {
      return true;
    }
    return false;
  }

  bool has_undirected_edge(int i, int j) {
    return g.at(i).at(j) && g.at(j).at(i);
  }
};

bool ci_test(const vector<vector<string>> &data, int node_x, int node_y, vector<int> Z, float ESS) {
  return true;
}

PDAG recursive_search(const vector<vector<string>> &data, PDAG &Gall, vector<int> Gs, vector<int> Gex, int N, float ESS) {
  int n_node = data.at(0).size();
  //stage 0: exit condition
  bool exitcondition = true;
  for (int i = 0; i < Gs.size(); i++) {
    if (Gall.predecessors(Gs.at(i)).size() > N){
      exitcondition = false;
      break;
    }
  }
  if (exitcondition) {
    return Gall;
  }
  //stage A1: Do CI tests for nodes between Gex and Gs and remove edge
  if (!(Gex.empty())){
    vector<vector<bool>>deletededges(n_node, vector<bool>(n_node, false));
    for (auto& node_y : Gs) {
      for (auto& node_x : Gex) {
        if (Gall.has_edge(node_x, node_y) && !(deletededges.at(node_x).at(node_y) || deletededges.at(node_y).at(node_x))) {
          vector<int> Z = Gall.predecessors(node_y);
          for (int i = 0; i < Z.size(); i++) {
            if (Z.at(i) == node_x) {
              Z.erase(Z.begin() + i);
            }
          } //erase node_x from Z
          if (Z.size() >= N) {
            foreach_comb(Z.size(), N, [&Z, &Gall, &data, &node_x, &node_y, &N, &ESS, &deletededges](int *indexes) { //error?
              //cout << indexes[0] << ',' << indexes[1] << endl;
              vector<int> selected_z;
              for (int j = 0; j < N; j++) {
                selected_z.push_back(Z.at(indexes[j]));
              }
              if (ci_test(data, node_x, node_y, selected_z, ESS)) {
                Gall.remove_edge(node_x, node_y);
                deletededges.at(node_x).at(node_y) = true;
                deletededges.at(node_y).at(node_x) = true;
                //transive_cut
              }
            });
          }
        }
      }
    }
  }
  

  return Gall;
}


PDAG RAI(const vector<vector<string>> &data, float ESS) {
  int n_node = data.at(0).size(); // number of nodes
  // initialize Gall, Gs, Gex
  PDAG Gall;
  vector<vector<bool>> gall(n_node, vector<bool>(n_node, true));
  for  (int i = 0; i < n_node; i++) {
          gall.at(i).at(i) = false;
      }
  Gall.g = gall; //Gall is complete graph
  vector<int> Gs(n_node, 1);// Gs contains all nodes 0 ~ n_node - 1
  vector<int> Gex;// Gex is empty

  PDAG Gend;
  Gend = recursive_search(data, Gall, Gs, Gex, 0, ESS);
  return Gend;
}


int main() {
  // dummy data,ESS
  vector<vector<string>> data = {{"a", "b", "c"}, {"a", "b", "c"}, {"a", "b", "c"}};
  float ESS = 1.0;
  
  PDAG Gend;
  Gend = RAI(data, ESS);
}
