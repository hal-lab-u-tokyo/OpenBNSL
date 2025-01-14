#define _GLIBCXX_DEBUG
using namespace std;
#include "structure_learning/PC.h"

#include <cassert>
#include <set>
#include <stack>
#include <string>
#include <vector>
// the following includes are for permutation and combination algorithms
#include <algorithm>
#include <functional>
#include <iostream>

// for gamma function
#include <cmath>

// input: data: np.ndarray,  shape: (n: number of variables, d: number of
// samples) output: leard PDAG
// Gall.at(i).at(j)==1 means there is an edge i -> j

template <typename T>
bool next_combination(const T first, const T last,
                      int k) {  // use sorted vector
  const T subset = first + k;
  // empty container | k = 0 | k == n
  if (first == last || first == subset || last == subset) {
    return false;
  }
  T src = subset;
  while (first != src) {
    src--;
    if (*src < *(last - 1)) {
      T dest = subset;
      while (*src >= *dest) {
        dest++;
      }
      iter_swap(src, dest);
      rotate(src + 1, dest + 1, last);
      rotate(subset, subset + (last - dest) - 1, last);
      return true;
    }
  }
  // restore
  rotate(first, subset, last);
  return false;
}

struct PDAG {
  vector<vector<bool>> g;
  // コンストラクタ
  PDAG() {
    // cout << "normal constructor called" << endl;
  }
  // コピーコンストラクタ
  PDAG(const PDAG &old) {
    // cout << "copy constructor called" << endl;
    g = old.g;
  }
  // 代入演算子
  PDAG &operator=(const PDAG &a) {
    if (this != &a) g = a.g;
    return *this;
  }
  // デストラクタ
  ~PDAG() = default;

  vector<int> successors(int i) {
    // return the list of successors of node i (include undirected edge)
    vector<int> succ;
    for (int j = 0; j < (int)g.size(); j++) {
      if (g.at(i).at(j)) {
        succ.push_back(j);
      }
    }
    return succ;
  }

  vector<int> predecessors(int i) {
    // return the list of predecessors of node i (include undirected edge)
    vector<int> pred;
    for (int j = 0; j < (int)g.size(); j++) {
      if (g.at(j).at(i)) {
        pred.push_back(j);
      }
    }
    return pred;
  }

  vector<int> neighbors(int i) {
    // return the list of neighbors {j} of node i (j -> i or i -> j)
    vector<int> neigh;
    for (int j = 0; j < (int)g.size(); j++) {
      if (g.at(j).at(i) || g.at(i).at(j)) {
        neigh.push_back(j);
      }
    }
    return neigh;
  }

  vector<int> undirected_neighbors(int i) {
    // return the list of undirected neighbors of node i
    vector<int> neigh;
    for (int j = 0; j < (int)g.size(); j++) {
      if (g.at(i).at(j) && g.at(j).at(i)) {
        neigh.push_back(j);
      }
    }
    return neigh;
  }

  void remove_edge(int i, int j) {
    // remove the edge i -> j
    g.at(i).at(j) = false;
  }

  void remove_edge_completedly(int i, int j) {
    // remove edge between i and j
    g.at(i).at(j) = false;
    g.at(j).at(i) = false;
  }

  void add_edge(int i, int j) { g.at(i).at(j) = true; }

  bool has_edge(int i, int j) { return g.at(i).at(j); }

  bool has_directed_edge(int i, int j) {
    if (g.at(i).at(j) && !g.at(j).at(i)) {
      return true;
    }
    return false;
  }

  bool has_undirected_edge(int i, int j) {
    return g.at(i).at(j) && g.at(j).at(i);
  }

  bool has_directed_path(
      int X,
      int Y) {  // check if there is a directed path from X to Y using DFS
    vector<int> visited(g.size(), 0);
    vector<int> stack;
    stack.push_back(X);
    while (!stack.empty()) {
      int node = stack.back();
      visited.at(node) = 1;
      stack.pop_back();
      if (node == Y) {
        return true;
      }
      for (auto &succ : successors(node)) {
        if (visited.at(succ) == 0 && has_directed_edge(node, succ)) {
          stack.push_back(succ);
        }
      }
    }
    return false;
  }
};

vector<vector<int>> state_count_PC(const vector<vector<int>> &data, int &node_x,
                                   vector<int> &parents,
                                   vector<int> &n_states) {
  if (parents.empty()) {
    int x = n_states.at(node_x);
    vector<vector<int>> counts(x, vector<int>(1, 0));
#pragma omp parallel
    {
      vector<int> temp(x, 0);
#pragma omp for
      for (int i = 0; i < (int)data.size(); i++) {
        temp.at(data.at(i).at(node_x)) += 1;
      }
      for (int j = 0; j < x; j++) {
#pragma omp atomic
        counts.at(j).at(0) += temp.at(j);
      }
    }
    return counts;
  } else {
    // return the state counts of X, Y | Z shape: countmap[state of child][state
    // of parents]
    int x = n_states.at(node_x);
    int y = 1;
    for (int i = 0; i < (int)parents.size(); i++) {
      y = y * n_states.at(parents.at(i));
    }
    vector<vector<int>> counts(x, vector<int>(y, 0));
    // count the number of each state
    int yy;
#pragma omp parallel private(yy)
    {
      vector<vector<int>> temp(x, vector<int>(y, 0));
#pragma omp for
      for (int i = 0; i < (int)data.size(); i++) {
        yy = data.at(i).at(parents.at(0));
        for (int j = 0; j < (int)parents.size(); j++) {
          yy = n_states.at(parents.at(j)) * yy + data.at(i).at(parents.at(j));
        }
        temp.at(data.at(i).at(node_x)).at(yy) += 1;
      }
      for (int j = 0; j < x; j++) {
        for (int k = 0; k < y; k++) {
#pragma omp atomic
          counts.at(j).at(k) += temp.at(j).at(k);
        }
      }
    }
    return counts;
  }
}

float localBDeuscore_PC(const vector<vector<int>> &data, int &node_x,
                        vector<int> &parents, double &ESS,
                        vector<int> &n_states) {
  // return log of the BDeu score of X, Y | Z
  double score = 0.0;
  if (parents.empty()) {
    // no parents
    vector<vector<int>> count;
    count = state_count_PC(data, node_x, parents, n_states);
    int r = count.size();  // number of states of node_x
    int n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count.at(k).at(0);
    }
    for (int k = 0; k < r; k++) {  // for each state of node_x
      score += lgamma(ESS / r + count.at(k).at(0)) - lgamma(ESS / r);
    }
    score += lgamma(ESS) - lgamma(ESS + n_i);
  } else {
    // have parents
    vector<vector<int>> count;
    count = state_count_PC(data, node_x, parents, n_states);
    int q = count.at(0).size();  // number of states of parents
    int r = count.size();        // number of states of node_x
    vector<float> n_ij(q, 0.0);
    for (int k = 0; k < r; k++) {
      for (int j = 0; j < q; j++) {
        n_ij.at(j) += count.at(k).at(j);
      }
    }
    for (int j = 0; j < q; j++) {    // for each state of parents
      for (int k = 0; k < r; k++) {  // for each state of node_x
        score +=
            lgamma(ESS / (r * q) + count.at(k).at(j)) - lgamma(ESS / (r * q));
        if (isnan(score)) {
          cout << ESS / (r * q) + count.at(k).at(j) << "score is nan" << endl;
        }
      }
      score += lgamma(ESS / q) - lgamma(ESS / q + n_ij.at(j));
    }
  }
  // calculate the score
  return score;
}

bool ci_test_PC(const vector<vector<int>> &data, int &node_x, int &node_y,
                vector<int> &Z, double &ESS, vector<int> &n_states) {
  // CI test for X _|_ Y | Z
  float independent_score = 0;
  float dependent_score = 0;
  float temp = localBDeuscore_PC(data, node_x, Z, ESS, n_states);
  independent_score += temp;
  independent_score += localBDeuscore_PC(data, node_y, Z, ESS, n_states);
  dependent_score += temp;
  vector<int> zplusx = Z;
  zplusx.push_back(node_x);
  dependent_score += localBDeuscore_PC(data, node_y, zplusx, ESS, n_states);
  if (independent_score > dependent_score) {
    // cout<< "CI independent:" <<node_x<<" _|_"<<node_y<<" |
    // "<<independent_score<<">"<<dependent_score<< endl;
    return true;
  } else {
    // cout<< "CI dependent:" <<node_x<<" _|_"<<node_y<<" |
    // "<<independent_score<<"<"<<dependent_score<< endl;
    return false;
  }
}

void orientation(PDAG &Gall, const vector<vector<int>> &data, double &ESS,
                 vector<int> &n_states) {
  /*
      orient edges in a PDAG to a maximally oriented graph.
      orient rules are based on rule 1~3 from Meek,C.:Causal Inference and
     Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in
     Artificial Inteligence (UAl-95),p.403-410 (195)
  */
  // for each X-Z-Y (X and Y is not adjecent), find V-structure and orient as X
  // -> Z <- Y
  for (int X = 0; X < (int)Gall.g.size(); X++) {
    for (auto &Z : Gall.undirected_neighbors(X)) {
      for (auto &Y : Gall.undirected_neighbors(Z)) {
        if (X != Y && !Gall.has_edge(X, Y) && !Gall.has_edge(Y, X)) {
          vector<int> z = {Z};
          // cout<< "V-structure think:" <<X<<"->"<<Z<<"<-"<<Y<< endl;
          if (!ci_test_PC(data, X, Y, z, ESS, n_states)) {
            Gall.remove_edge(Z, X);
            Gall.remove_edge(Z, Y);
            // cout<< "V-structure found:" <<X<<"->"<<Z<<"<-"<<Y<< endl;
          }
        }
      }
    }
  }
  bool flag = true;
  while (flag) {
    flag = false;
    // Rule 1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
    for (int X = 0; X < (int)Gall.g.size(); X++) {
      for (auto &Y : Gall.successors(X)) {
        if (!Gall.has_edge(Y, X) && Gall.has_edge(X, Y)) {
          for (auto &Z : Gall.undirected_neighbors(Y)) {
            if (!Gall.has_edge(X, Z) && !Gall.has_edge(Z, X) && Z != X) {
              Gall.remove_edge(Z, Y);
              // cout<< "R1:" <<Y<<"->"<<Z<< endl;
              flag = true;
            }
          }
        }
      }
    }
    // Rule 2: X - Y and if there is a directed path from X to Y, then X -> Y
    for (int X = 0; X < (int)Gall.g.size(); X++) {
      for (auto &Y : Gall.undirected_neighbors(X)) {
        if (Gall.has_directed_path(X, Y)) {
          Gall.remove_edge(Y, X);
          // cout<< "R2:" <<X<<"->"<<Y<< endl;
          flag = true;
        }
      }
    }
    // Rule 3: for each X->W<-Z X-Y-Z Y-W, orient Y->W
    for (int X = 0; X < (int)Gall.g.size(); X++) {
      for (auto &Y : Gall.undirected_neighbors(X)) {
        for (auto &Z : Gall.undirected_neighbors(Y)) {
          if (Z != X) {  // X-Y-Z
            for (auto &W : Gall.undirected_neighbors(Y)) {
              if (W != X && W != Z && Gall.has_directed_edge(X, W) &&
                  Gall.has_directed_edge(Z, W)) {
                Gall.remove_edge(W, Y);
                // cout<< "R3:" <<Y<<"->"<<W<< endl;
                flag = true;
              }
            }
          }
        }
      }
    }
  }

  return;
}

PDAG PCsearch(const vector<vector<int>> &data, PDAG &Gall, double &ESS,
              vector<int> &n_states) {
  int n_node = data.at(0).size();

  // stage 1: Do CI tests between nodes and remove edge (undirected graph)
  int t = 0;
  while (t < n_node - 2) {
    for (int node_x = 1; node_x < n_node; node_x++) {
      for (int node_y = 0; node_y < node_x; node_y++) {
        if (t == 0) {
          vector<int> S;
          if (ci_test_PC(data, node_x, node_y, S, ESS, n_states)) {
            Gall.remove_edge(node_x, node_y);
            Gall.remove_edge(node_y, node_x);
          }
        } else {
          vector<int> S;
          for (auto &node_z : Gall.predecessors(node_x)) {
            if (Gall.has_edge(node_y, node_z)) {
              S.push_back(node_z);
            }
          }
          if ((int)S.size() >= t) {
            do {
              vector<int> selected_z;
              for (int j = 0; j < t; j++) {
                selected_z.push_back(S.at(j));
              }
              if (ci_test_PC(data, node_x, node_y, selected_z, ESS, n_states)) {
                Gall.remove_edge(node_x, node_y);
                Gall.remove_edge(node_y, node_x);
              }
            } while (next_combination(S.begin(), S.end(), t));
          }
        }
      }
    }
    t++;
  }
  // stage 2: orient edges
  orientation(Gall, data, ESS, n_states);

  return Gall;
}

// void RAI(const vector<vector<int>> data, double ESS,vector<vector<bool>> g) {
py::array_t<bool> PC(py::array_t<int> data, py::array_t<int> n_states,
                     double ESS) {
  // translate imput data to c++ vector(this is not optimal but I don't know how
  // to use pybind11::array_t)
  py::buffer_info buf_data = data.request(), buf_states = n_states.request();
  const int *__restrict__ prt_data = static_cast<int *>(buf_data.ptr);
  const int *__restrict__ prt_states = static_cast<int *>(buf_states.ptr);
  size_t n_data = buf_data.shape[0],
         n_node = buf_data.shape[1];  // number of nodes
  vector<vector<int>> data_vec(n_data, vector<int>(n_node));
  vector<int> n_states_vec(n_node, 0);
  for (size_t i = 0; i < n_data; i++) {
    for (size_t j = 0; j < n_node; j++) {
      data_vec.at(i).at(j) = prt_data[i * n_node + j];
    }
  }
  for (size_t i = 0; i < n_node; i++) {
    n_states_vec.at(i) = prt_states[i];
  }
  auto endg = py::array_t<bool>({n_node, n_node});
  py::buffer_info buf_endg = endg.request();
  bool *__restrict__ prt_endg = static_cast<bool *>(buf_endg.ptr);

  // initialize Gall (complete undirected graph)
  PDAG Gall;
  vector<vector<bool>> gall(n_node, vector<bool>(n_node, true));
  for (int i = 0; i < (int)n_node; i++) {
    gall.at(i).at(i) = false;
  }
  Gall.g = gall;  // Gall is complete graph
  PDAG Gend;
  // vector<vector<string>> state_list = get_state_list(data);
  Gend = PCsearch(data_vec, Gall, ESS, n_states_vec);

  // translate Gend to py::array_t (this is not optimal but I don't know how to
  // use pybind11::array_t)
  for (size_t i = 0; i < n_node; i++) {
    for (size_t j = 0; j < n_node; j++) {
      prt_endg[i * n_node + j] = Gend.g.at(i).at(j);
    }
  }
  return endg;
}
