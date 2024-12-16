#define _GLIBCXX_DEBUG
using namespace std;
#include "structure_learning/RAI.h"
#include <string>
#include <set>
#include<vector>
#include<stack>
#include<cassert>
//the following includes are for permutation and combination algorithms
#include <algorithm>
#include <functional>
#include <iostream>

//for gamma function
#include <cmath>



// input: data: np.ndarray,  shape: (n: number of variables, d: number of
// samples) output: leard PDAG
// Gall.at(i).at(j)==1 means there is an edge i -> j
// c++ 17 を仮定]


void dfs(int now, vector<vector<int>> &G, vector<bool> &visited, vector<int> &order) {
    //cout<<now<<endl;
    visited[now] = true;
    for (int i = 0; i < (int)G.at(now).size(); i++) {
        if (!visited[G.at(now).at(i)]) {
            dfs(G.at(now).at(i), G, visited, order);
        }
    }
    order.push_back(now);
}

void rdfs(int now, vector<vector<int>> &rG, vector<bool> &visited, int k, vector<int> &comp) {
    //cout<<now<<endl;
    visited[now] = true;
    comp[now] = k;
    for (int i = 0; i < (int)rG.at(now).size(); i++) {
        if (!visited[rG.at(now).at(i)]) {
            rdfs(rG.at(now).at(i), rG, visited, k, comp);
        }
    }
}

vector<int> SCC(vector<vector<int>> &G, vector<vector<int>> &rG){
    int n_node = G.size();
    vector<bool> visited(n_node, false);
    vector<int> comp(n_node, 0);
    vector<int> order;
    order.clear();
    //dfs
    for (int i = 0; i < n_node; i++) {
        if (!visited[i]) {
            dfs(i, G, visited, order);
        }
    }
    fill(visited.begin(), visited.end(), false);
    //dfs2
    int k = 0;
    for (int i = n_node - 1; i >= 0; i--) {
        if (!visited[order[i]]) {
            rdfs(order.at(i), rG, visited, k, comp);
            k++;
        }
    }
    return comp;
}  //lowest topological order -> highest number in comp

vector<vector<int>> adjmat2listmat(vector<vector<bool>> &adjmat){//隣接行列表現2隣接リスト表現
    int n_node = adjmat.size();
    vector<vector<int>> listmat(n_node);
    for (int i = 0; i < n_node; i++) {
        for (int j = 0; j < n_node; j++) {
            if (adjmat.at(i).at(j)) {
                listmat.at(i).push_back(j); //i to j
            }
        }
    }
    return listmat;
}

vector<vector<int>> adjmat2listmat_reverse(vector<vector<bool>> &adjmat){//隣接行列表現2隣接リスト表現(reversed)
    int n_node = adjmat.size();
    vector<vector<int>> listmat(n_node);
    for (int i = 0; i < n_node; i++) {
        for (int j = 0; j < n_node; j++) {
            if (adjmat.at(i).at(j)) {
                listmat.at(j).push_back(i); //reverse edge j to i
            }
        }
    }
    return listmat;
}

vector<vector<bool>> make_gs_graph(vector<vector<bool>> &Gallgraph, vector<int> &gs){
    vector<vector<bool>> gsgraph(gs.size(), vector<bool>(gs.size(), false));
    for (int i = 0; i < gs.size(); i++){
        for (int j = 0; j < gs.size(); j++){
            gsgraph.at(i).at(j) = Gallgraph.at(gs.at(i)).at(gs.at(j));
        }
    }
    return gsgraph;
}

template <typename T> bool next_combination(const T first, const T last, int k) { //use sorted vector
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
    //cout << "normal constructor called" << endl;
  }
  // コピーコンストラクタ
  PDAG(const PDAG &old) {
    //cout << "copy constructor called" << endl;
    g = old.g;
  }
  // 代入演算子
  PDAG& operator=(const PDAG& a) {
    if (this != &a) g = a.g;
    return *this;
  }
  // デストラクタ
  ~PDAG() = default;

  vector<int> successors(int i) {
    // return the list of successors of node i (include undirected edge)
    vector<int> succ;
    for (int j = 0; j < g.size(); j++) {
      if (g.at(i).at(j)) {
        succ.push_back(j);
      }
    }
    return succ;
  }

  vector<int> predecessors(int i) {
    // return the list of predecessors of node i (include undirected edge)
    vector<int> pred;
    for (int j = 0; j < g.size(); j++) {
      if (g.at(j).at(i)) {
        pred.push_back(j);
      }
    }
    return pred;
  }

  vector<int> neighbors(int i) { 
    //return the list of neighbors {j} of node i (j -> i or i -> j)
    vector<int> neigh;
    for (int j = 0; j < g.size(); j++) {
      if (g.at(j).at(i) || g.at(i).at(j)) {
        neigh.push_back(j);
      }
    }
    return neigh;
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

  void remove_edge(int i, int j) {
    // remove the edge i -> j
    g.at(i).at(j) = false;
  }

  void remove_edge_completedly(int i, int j) {
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

  bool has_directed_path(int X, int Y){ //check if there is a directed path from X to Y using DFS
    vector<int> visited(g.size(), 0);
    vector<int> stack;
    stack.push_back(X);
    while (!stack.empty()){
      int node = stack.back();
      visited.at(node) = 1;
      stack.pop_back();
      if (node == Y){
        return true;
      }
      for (auto& succ : successors(node)){
        if (visited.at(succ) == 0 && has_directed_edge(node, succ)){
          stack.push_back(succ);
        }
      }
    }
    return false;
  }
};

void  order_grouping(PDAG &Gall, vector<int> &Gs, vector<int> &Gd, vector<vector<int>> &g_subs) {
  vector<vector<bool>> adjmat = make_gs_graph(Gall.g, Gs);
  vector<vector<int>> listmat = adjmat2listmat(adjmat);
  vector<vector<int>> rlistmat = adjmat2listmat_reverse(adjmat);
  vector<int> comp = SCC(listmat, rlistmat);
  bool flag = true;
  int order = 0;
  while(flag){
    flag = false;
    for (int i = 0; i < comp.size(); i++){
      if (comp.at(i) == order){ // node Gs.at(i) is in the order-th group
        g_subs.push_back(vector<int>());
        g_subs.at(order).push_back(Gs.at(i)); // node Gs.at(i) is in the order-th group(g_subs.at(order))
        flag = true;
      }
    }
    order = order + 1;
  }
  Gd = g_subs.back();
  g_subs.pop_back();
  return;
}

vector<vector<int>> state_count(const vector<vector<int>> &data, int &node_x, vector<int> &parents, vector<int> &n_states) {
  if(parents.empty()) {
    int x = n_states.at(node_x);
    vector<vector<int>> counts(x, vector<int>(1, 0));
    #pragma omp parallel for
      for(int i = 0; i < data.size(); i++) {
        counts.at(data.at(i).at(node_x)).at(0) += 1;
      }
    return counts;
  } else {
    //return the state counts of X, Y | Z shape: countmap[state of child][state of parents]
    int x = n_states.at(node_x);
    int y = 1;
    for (int i = 0; i < parents.size(); i++) {
      y = y * n_states.at(parents.at(i));
    }
    vector<vector<int>> counts(x, vector<int>(y, 0));
    //count the number of each state
    #pragma omp parallel for
      for(int i = 0; i < data.size(); i++) {
        int yy;
        for (int j = 0; j < parents.size(); j++) {
          if (j == 0) {
            yy = data.at(i).at(parents.at(j));
          }
          else {
            yy = n_states.at(parents.at(j)) * yy + data.at(i).at(parents.at(j));
          }
        }
        counts.at(data.at(i).at(node_x)).at(yy) += 1;
      }
    return counts;
  }
}

float localBDeuscore(const vector<vector<int>> &data, int &node_x, vector<int> &parents, double &ESS, vector<int> &n_states) {
  //return log of the BDeu score of X, Y | Z
  double score = 0.0;
  if (parents.empty()) {
    //no parents
    vector<vector<int>> count;
    count = state_count(data, node_x, parents, n_states);
    int r = count.size(); //number of states of node_x
    int n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count.at(k).at(0);
    }
    for (int k = 0; k < r; k++) { //for each state of node_x
        score += lgamma(ESS / r + count.at(k).at(0)) - lgamma(ESS/r);
    }
    score += lgamma(ESS) - lgamma(ESS + n_i);
  }
  else {
    //have parents
    vector<vector<int>> count;
    count = state_count(data, node_x, parents, n_states);
    int q = count.at(0).size(); //number of states of parents
    int r = count.size(); //number of states of node_x
    vector<float> n_ij(q, 0.0);
    for (int k = 0; k < r; k++) {
      for(int j = 0; j < q; j++) {
        n_ij.at(j) += count.at(k).at(j);
      }
    }
    for (int j = 0; j < q; j++) { //for each state of parents
      for (int k = 0; k < r; k++) { //for each state of node_x
        score += lgamma(ESS / (r * q) + count.at(k).at(j)) - lgamma(ESS/(r * q));
        if (isnan(score)) {
          cout <<ESS / (r * q) + count.at(k).at(j)<< "score is nan" << endl;
        }
      }
      score += lgamma(ESS / q) - lgamma(ESS / q + n_ij.at(j));
    }
  }
    //calculate the score
  return score;
}

bool ci_test(const vector<vector<int>> &data, int &node_x, int &node_y, vector<int> &Z, double &ESS, vector<int> &n_states) {
  //CI test for X _|_ Y | Z
  float independent_score = 0;
  float dependent_score = 0;
  float temp = localBDeuscore(data, node_x, Z, ESS, n_states);
  independent_score += temp;
  independent_score += localBDeuscore(data, node_y, Z, ESS, n_states);
  dependent_score += temp;
  vector<int> zplusx = Z;
  zplusx.push_back(node_x);
  dependent_score += localBDeuscore(data, node_y, zplusx, ESS, n_states);
  if(independent_score > dependent_score){
    cout<< "CI independent:" <<node_x<<" _|_"<<node_y<<" | "<<independent_score<<">"<<dependent_score<< endl;
    return true;
  }else{
    cout<< "CI dependent:" <<node_x<<" _|_"<<node_y<<" | "<<independent_score<<"<"<<dependent_score<< endl;
    return false;
  }
}

// vector<vector<string>> get_state_list(const vector<vector<string>> &data){
//   vector<vector<string>> state_list(data.at(0).size());
//   for (int i = 0; i < data.at(0).size(); i++) { //node
//     for (int j = 0; j < data.size(); j++) { //n_sample
//       if(find(state_list.at(i).begin(), state_list.at(i).end(), data.at(j).at(i)) == state_list.at(i).end()){
//         state_list.at(i).push_back(data.at(j).at(i));
//       }
//     }
//   }
//   return state_list;
// }

// vector<vector<int>> data_string2int(const vector<vector<string>> &data) {
//   vector<vector<int>> data_int(data.size(), vector<int>(data.at(0).size()));
//   for (int i = 0; i < data.size(); i++) {
//     for (int j = 0; j < data.at(0).size(); j++) {
//       data_int.at(i).at(j) = find(state_list.at(j).begin(), state_list.at(j).end(), data.at(i).at(j)) - state_list.at(j).begin();
//     }
//   }
//   return data_int;

// }

void orientation_A2(PDAG &Gall, vector<int> &Gs, vector<vector<bool>> &deletededges) {
  /*
      orient edges in a PDAG to a maximally oriented graph.
      orient rules are based on rule 1~3 from Meek,C.:Causal Inference and Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in Artificial Inteligence (UAl-95),p.403-410 (195)
      in this stage (stage A2 in Yehezkel and Lerner(2009)), only rule 1 is applied because only X -> Y - Z shape is created in stage A1 (X-Z removed).
  */
  //Rule 1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
  for (int i = 0; i < Gall.g.size(); i++){
    for (int j = 0; j < Gall.g.size(); j++){
      if (deletededges.at(i).at(j) || deletededges.at(j).at(i)){
        int X = i;
        int Z = j;
        for (auto& Y : Gall.undirected_neighbors(Z)) {
          if (Gall.has_directed_edge(X, Y)) {
            Gall.remove_edge(Z, Y);
            cout<< "A2_removed:" <<Y<<"->"<<Z<< endl;
          }
        }
      }
    }
  }
  return;
}

void orientation_B2(PDAG &Gall, vector<int> &Gs, vector<vector<bool>> &deletededges, const vector<vector<int>> &data, double &ESS, vector<int> &n_states) {
  /*
      orient edges in a PDAG to a maximally oriented graph.
      orient rules are based on rule 1~3 from Meek,C.:Causal Inference and Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in Artificial Inteligence (UAl-95),p.403-410 (195)
  */
  //for each X-Z-Y (X and Y is not adjecent), find V-structure and orient as X -> Z <- Y
  for (auto& X : Gs) {
    for (auto& Z : Gall.undirected_neighbors(X)) {
      for (auto& Y : Gall.undirected_neighbors(Z)) {
        if (X != Y && !Gall.has_edge(X, Y) && !Gall.has_edge(Y, X)) {
          vector<int> z = {Z};
          cout<< "V-structure think:" <<X<<"->"<<Z<<"<-"<<Y<< endl;
          if (!ci_test(data, X, Y, z, ESS, n_states)){
            Gall.remove_edge(Z, X);
            Gall.remove_edge(Z, Y);
            cout<< "V-structure found:" <<X<<"->"<<Z<<"<-"<<Y<< endl;
          }
        }
      }
    }
  }
  bool flag = true;
  while (flag){
    flag = false;
    //Rule 1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
    for (auto& X : Gs) {
      for (auto& Y : Gall.successors(X)) {
        if (!Gall.has_directed_edge(Y, X)) {
          for (auto& Z : Gall.undirected_neighbors(Y)) {
            if (!Gall.has_edge(X, Z) && !Gall.has_edge(Z, X) && Z != X) {
              Gall.remove_edge(Z, Y);
              cout<< "R1:" <<Y<<"->"<<Z<< endl;
              flag = true;
            }
          }
        }
      }
    }
    //Rule 2: X - Y and if there is a directed path from X to Y, then X -> Y
    for (auto& X : Gs) {
      for (auto& Y : Gall.undirected_neighbors(X)) {
        if (Gall.has_directed_path(X, Y)) {
          Gall.remove_edge(Y, X);
          cout<< "R2:" <<X<<"->"<<Y<< endl;
          flag = true;
        }
      }
    }
    //Rule 3: for each X->W<-Z X-Y-Z Y-W, orient Y->W
    for (auto &X : Gs){
      for (auto &Y : Gall.undirected_neighbors(X)){
        for (auto &Z : Gall.undirected_neighbors(Y)){
          if (Z != X){ //X-Y-Z
            for (auto &W : Gall.undirected_neighbors(Y)){
              if (W != X && W != Z && Gall.has_directed_edge(X, W) && Gall.has_directed_edge(Z, W)){
                Gall.remove_edge(W, Y);
                cout<< "R3:" <<Y<<"->"<<W<< endl;
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

PDAG recursive_search(const vector<vector<int>> &data, PDAG &Gall, vector<int> Gs, vector<int> Gex, int N, double &ESS, vector<int> &n_states) {
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
  //stage A1: Do CI tests for nodes between Gex and Gs and remove edge   (bug?)
  vector<vector<bool>>deletededges(n_node, vector<bool>(n_node, false));
  if (!(Gex.empty())){
    for (auto& node_y : Gs) {
      for (auto& node_x : Gex) {
        if (Gall.has_edge(node_x, node_y) && !(deletededges.at(node_x).at(node_y) || deletededges.at(node_y).at(node_x))) {
          vector<int> Z = Gall.predecessors(node_y);
          for (int i = 0; i < Z.size(); i++) {
            if (Z.at(i) == node_x) {
              Z.erase(Z.begin() + i);
            }
          } //erase node_x from Z
          sort(Z.begin(), Z.end());
          if (Z.size() >= N) {
            do {
                  vector<int> selected_z;
                  for (int j = 0; j < N; j++) {
                    selected_z.push_back(Z.at(j));
                  }
                  if (ci_test(data, node_x, node_y, selected_z, ESS, n_states)) {
                    Gall.remove_edge(node_x, node_y);
                    deletededges.at(node_x).at(node_y) = true;
                    deletededges.at(node_y).at(node_x) = true;
                    cout<< "A1_removed:" <<node_x<<"-"<<node_y<< endl;
                    //transive_cut();
                  }
              } while(next_combination(Z.begin(), Z.end(), N));
          }
        }
      }
    }
  }
  //stage A2: orient edges in Gs using "smart" orientation rules R1
  orientation_A2(Gall, Gs, deletededges);

  //stage B1: Do CI tests for nodes between Gs and Gs and remove edge
  for (auto& node_y : Gs) {
    for (auto& node_x : Gall.neighbors(node_y)) {
      if (N == 0) {
        vector<int> S;
        if (ci_test(data, node_x, node_y, S, ESS, n_states)) {
          Gall.remove_edge(node_x, node_y);
          Gall.remove_edge(node_y, node_x);
          deletededges.at(node_x).at(node_y) = true;
          deletededges.at(node_y).at(node_x) = true;
          cout<< "B1_removed:" <<node_x<<"-"<<node_y<< endl;
          //transive_cut();
        }
      }else{
        vector<int> S = Gall.predecessors(node_y);
        auto newEnd = remove(S.begin(), S.end(), node_x);
        S.erase(newEnd, S.end()); // remove node_x from S  正しい？
        if (S.size() >= N) {
          do {
            vector<int> selected_z;
            for (int j = 0; j < N; j++) {
              selected_z.push_back(S.at(j));
            }
            if (ci_test(data, node_x, node_y, selected_z, ESS, n_states)) {
              Gall.remove_edge(node_x, node_y);
              Gall.remove_edge(node_y, node_x);
              deletededges.at(node_x).at(node_y) = true;
              deletededges.at(node_y).at(node_x) = true;
              cout<< "B1_removed:" <<node_x<<"-"<<node_y<< endl;
              //transive_cut();
            }
          } while(next_combination(S.begin(), S.end(), N));
        }
      }
    }
  }
  
  //stage B2: orient edges in Gs using orientation rules R1~R3
  orientation_B2(Gall, Gs, deletededges, data, ESS, n_states);
  //stage B3: Group the nodes having the lowest topological order into a descendant substructure Gd
  vector<int> Gd;
  vector<vector<int>> g_subs;
  order_grouping(Gall, Gs, Gd, g_subs);
  //stage C: Ancestorsub-structure decomposition
  vector<int> Gexd;
  for (auto& Gex_i : g_subs) {
    for (auto& node : Gex_i) {
      Gexd.push_back(node);
    }
    sort(Gexd.begin(), Gexd.end());
    Gall = recursive_search(data, Gall, Gex_i, Gex, N + 1, ESS, n_states);
  }
  //stage D:  Descendantsub-structure decomposition
  return recursive_search(data, Gall, Gd, Gexd, N + 1, ESS, n_states);
}

//void RAI(const vector<vector<int>> data, double ESS,vector<vector<bool>> g) {
py::array_t<bool> RAI(py::array_t<int> data, py::array_t<int> n_states, double ESS) {
  //translate imput data to c++ vector(this is not optimal but I don't know how to use pybind11::array_t) 
  py::buffer_info buf_data = data.request(), buf_states = n_states.request();
  const int* __restrict__ prt_data = static_cast<int*>(buf_data.ptr);
  const int* __restrict__ prt_states = static_cast<int*>(buf_states.ptr);
  size_t n_data = buf_data.shape[0], n_node = buf_data.shape[1]; // number of nodes
  vector<vector<int>> data_vec(n_data, vector<int>(n_node));
  vector<int> n_states_vec(n_node, 0);
  for(size_t i = 0; i< n_data; i++){
    for(size_t j = 0; j< n_node; j++){
      data_vec.at(i).at(j) = prt_data[i * n_node + j];
    }
  }
  for (size_t i = 0; i < n_node; i++) {
    n_states_vec.at(i) = prt_states[i];
  }
  auto endg = py::array_t<bool>({n_node, n_node});
  py::buffer_info buf_endg = endg.request();
  bool* __restrict__ prt_endg = static_cast<bool*>(buf_endg.ptr);
  
  // initialize Gall, Gs, Gex
  PDAG Gall;
  vector<vector<bool>> gall(n_node, vector<bool>(n_node, true));
  for  (int i = 0; i < n_node; i++) {
          gall.at(i).at(i) = false;
      }
  Gall.g = gall; //Gall is complete graph
  vector<int> Gs(n_node, -1);
  for (int i = 0; i < n_node; i++) {
    Gs.at(i) = i;
  }// Gs contains all nodes 0 ~ n_node - 1
  vector<int> Gex;// Gex is empty
  PDAG Gend;
  // vector<vector<string>> state_list = get_state_list(data);
  Gend = recursive_search(data_vec, Gall, Gs, Gex, 0, ESS, n_states_vec);

  //translate Gend to py::array_t (this is not optimal but I don't know how to use pybind11::array_t)
  for (size_t i = 0; i < n_node; i++) {
    for (size_t j = 0; j < n_node; j++) {
      prt_endg[i * n_node + j] = Gend.g.at(i).at(j);
    }
  }
  return endg;
}

