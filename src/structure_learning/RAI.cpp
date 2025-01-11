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
    //////////////cout<<now<<endl;
    visited[now] = true;
    for (int i = 0; i < (int)G.at(now).size(); i++) {
        if (!visited[G.at(now).at(i)]) {
            dfs(G.at(now).at(i), G, visited, order);
        }
    }
    order.push_back(now);
}

void rdfs(int now, vector<vector<int>> &rG, vector<bool> &visited, int k, vector<int> &comp) {
    //////////////cout<<now<<endl;
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
    //////////////cout << "normal constructor called" << endl;
  }
  // コピーコンストラクタ
  PDAG(const PDAG &old) {
    //////////////cout << "copy constructor called" << endl;
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

  bool has_path(int X, int Y){ //check if there is a directed path from X to Y using DFS
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
        if (visited.at(succ) == 0 && has_edge(node, succ)){
          stack.push_back(succ);
        }
      }
    }
    return false;
  }

   bool has_connection(int X, int Y){ //check if there is connection between X and Y using DFS
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
        if (visited.at(succ) == 0 && (has_edge(node, succ) || has_edge(succ, node))){
          stack.push_back(succ);
        }
      }
    }
    return false;
  } 
};

void  order_grouping(PDAG &Gall, vector<int> &Gs, vector<int> &Gd, vector<vector<int>> &g_ex_connection) {
  vector<vector<bool>> adjmat = make_gs_graph(Gall.g, Gs);
  vector<vector<int>> listmat = adjmat2listmat(adjmat);
  vector<vector<int>> rlistmat = adjmat2listmat_reverse(adjmat);
  vector<int> comp = SCC(listmat, rlistmat);
  PDAG Gss;
  Gss.g = make_gs_graph(Gall.g, Gs);
  bool flag = true;
  int order = 0;
  vector<vector<int>> g_subs;
  while(flag){
    flag = false;
    for (int i = 0; i < comp.size(); i++){
      if (comp.at(i) == order){ // node Gs.at(i) is in the order-th group
        if (flag == false){
          g_subs.push_back(vector<int>());
        }
        g_subs.at(order).push_back(i); // g_subsはGsのindexを格納
        flag = true;
      }
    }
    order = order + 1;
  }
  ////cout<< "comp:";
  for (int i = 0; i < comp.size(); i++){
    ////cout<< comp.at(i)<<", ";
  }
  ////cout<<endl;


  //子孫部分集合の分離
  //各クラスタについてその全てのノードが他のクラスタに子ノードを持たないもの(クラスタ内の一つのノードを持ってきたときに他のクラスタのどれか一つのノードへのpathがないもの)を判定
  vector<bool> is_Gd(g_subs.size(), true);
  for (int i = 0; i < g_subs.size(); i++){
    for (int j = 0; j < g_subs.size(); j++){
      if (i != j){
        if(Gss.has_path(g_subs.at(i).at(0), g_subs.at(j).at(0))){ //index gss と gall の不一致
          is_Gd.at(i) = false;
          break;
        }
      }
    }
  }
  //それら全てをG_dとする 
  vector<int> Gex;
  for (int i = 0; i < g_subs.size(); i++){
    if (is_Gd.at(i)){
      for (int j = 0; j < g_subs.at(i).size(); j++){
        Gd.push_back(g_subs.at(i).at(j));
      }
    }else{
      for (int j = 0; j < g_subs.at(i).size(); j++){
        Gex.push_back(g_subs.at(i).at(j));
      }
    }
  }
  //G_dが全てなら終了
  bool temp = true;
  for (int i = 0; i < is_Gd.size(); i++){
    if (!is_Gd.at(i)){
      temp = false;
      break;
    }
  }
  if (temp){
    return;
  }

  //G_d以外について，クラスタ0から連結しているクラスタをまとめていく
  PDAG Gss_ex;
  Gss_ex.g = make_gs_graph(Gss.g, Gex);
  for (int i = 0; i < Gss_ex.g.size(); i++){
    for (int j = 0; j < Gss_ex.g.size(); j++){
      if (Gss_ex.g.at(i).at(j)){
        Gss_ex.g.at(j).at(i) = true;
      }
    }
  } //adjmat_ex is undirected graph
  vector<int> Gex_reverseindexlist(Gss.g.size(), -1); //.at(Gs_node_index) = Gex_node_index
  for (int i = 0; i < Gex.size(); i++){
    Gex_reverseindexlist.at(Gex.at(i)) = i;
  }

  vector<bool> flag_ex_vec(g_subs.size(), true);
  for (int i = 0; i < g_subs.size(); i++){
    if(flag_ex_vec.at(i) && !is_Gd.at(i)){
      flag_ex_vec.at(i) = false;
      g_ex_connection.push_back(vector<int>());
      for (int l = 0; l < g_subs.at(i).size(); l++){
        g_ex_connection.back().push_back(g_subs.at(i).at(l));
      }
      for (int j = 0; j < g_subs.size(); j++){
        if (i != j && !is_Gd.at(j) && !g_subs.at(i).empty() && !g_subs.at(j).empty()){
          if (Gss_ex.has_connection(Gex_reverseindexlist.at(g_subs.at(i).at(0)), Gex_reverseindexlist.at(g_subs.at(j).at(0)))){
            for (int k = 0; k < g_subs.at(j).size(); k++){
              g_ex_connection.back().push_back(g_subs.at(j).at(k));
            }
            flag_ex_vec.at(j) = false;
          }
        }
      }
      sort(g_ex_connection.back().begin(), g_ex_connection.back().end());
    }
  }
  //Gs, g_ex_connectionのindexをGallのindexに変換
  for (int i = 0; i < Gd.size(); i++){
    Gd.at(i) = Gs.at(Gd.at(i));
  }
  for (int i = 0; i < g_ex_connection.size(); i++){
    for (int j = 0; j < g_ex_connection.at(i).size(); j++){
      g_ex_connection.at(i).at(j) = Gs.at(g_ex_connection.at(i).at(j));
    }
  }
  return;
}

vector<vector<int>> state_count(const vector<vector<int>> &data, vector<int> &children, vector<int> &parents, vector<int> &n_states) {
  if (children.size() == 1){
    int node_x = children.at(0);
    if(parents.empty()) {
      int x = n_states.at(node_x);
      vector<vector<int>> counts(x, vector<int>(1, 0));
      #pragma omp parallel
        {
          vector<int> temp(x, 0);
          #pragma omp for
            for(int i = 0; i < data.size(); i++) {
              temp.at(data.at(i).at(node_x)) += 1;
            }
          for(int j = 0; j < x; j++){
            #pragma omp atomic
            counts.at(j).at(0) += temp.at(j);
          }
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
      int yy;
      #pragma omp parallel private(yy) 
        {
          vector<vector<int>> temp(x, vector<int>(y, 0));
          #pragma omp for
            for(int i = 0; i < data.size(); i++) {
              for (int j = 0; j < parents.size(); j++) {
                if (j == 0) {
                  yy = data.at(i).at(parents.at(j));
                }
                else {
                  yy = n_states.at(parents.at(j)) * yy + data.at(i).at(parents.at(j));
                }
              }
              temp.at(data.at(i).at(node_x)).at(yy) += 1;
            }
          for(int j = 0; j < x; j++){
            for(int k = 0; k < y; k++){
              #pragma omp atomic
              counts.at(j).at(k) += temp.at(j).at(k);
            }
          }
        }
      return counts;
    }
  }else{
    if(parents.empty()) {
      int xc = n_states.at(children.at(0));
      int yc = n_states.at(children.at(1));
      int len;
      len = xc * yc;
      vector<vector<int>> counts(len, vector<int>(1, 0));
      #pragma omp parallel
        {
          vector<int> temp(len, 0);
          #pragma omp for
            for(int i = 0; i < data.size(); i++) {
              temp.at(data.at(i).at(children.at(0)) * yc + data.at(i).at(children.at(1))) += 1;
            }
          for(int j = 0; j < len; j++){
            #pragma omp atomic
            counts.at(j).at(0) += temp.at(j);
          }
        }
      return counts;
    } else {
      //return the state counts of X, Y | Z shape: countmap[state of child][state of parents]
      int xc = n_states.at(children.at(0));
      int yc = n_states.at(children.at(1));
      int len;
      len = xc * yc;
      int y = 1;
      for (int i = 0; i < parents.size(); i++) {
        y = y * n_states.at(parents.at(i));
      }
      vector<vector<int>> counts(len, vector<int>(y, 0));
      //count the number of each state
      int yy;
      #pragma omp parallel private(yy) 
        {
          vector<vector<int>> temp(len, vector<int>(y, 0));
          #pragma omp for
            for(int i = 0; i < data.size(); i++) {
              for (int j = 0; j < parents.size(); j++) {
                if (j == 0) {
                  yy = data.at(i).at(parents.at(j));
                }
                else {
                  yy = n_states.at(parents.at(j)) * yy + data.at(i).at(parents.at(j));
                }
              }
              temp.at(data.at(i).at(children.at(0)) * yc + data.at(i).at(children.at(1))).at(yy) += 1;
            }
          for(int j = 0; j < len; j++){
            for(int k = 0; k < y; k++){
              #pragma omp atomic
              counts.at(j).at(k) += temp.at(j).at(k);
            }
          }
        }
      return counts;
    }
  }
}

float natori_independent_score(const vector<vector<int>> &data, int &node_x, int &node_y, vector<int> &parents, vector<int> &n_states, float &ESS) {
  //return log of the BDeu score of X, Y | Z
  double score = 0.0;
  double alpha = 0.5; // hyperparameter
  vector<int> node_x_vec(1, node_x);
  vector<int> node_y_vec(1, node_y);
  if (parents.empty()) {
    //no parents
    vector<vector<int>> count;
    count = state_count(data, node_x_vec, parents, n_states);
    int r = count.size(); //number of states of node_x
    int n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count.at(k).at(0);
    }
    for (int k = 0; k < r; k++) { //for each state of node_x
        score += lgamma(alpha + count.at(k).at(0)) - lgamma(alpha);
    }
    score += lgamma(r * alpha) - lgamma(r * alpha + n_i);
    vector<vector<int>> count2;
    count2 = state_count(data, node_y_vec, parents, n_states);
    r = count2.size(); //number of states of node_x
    n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count2.at(k).at(0);
    }
    for (int k = 0; k < r; k++) { //for each state of node_x
        score += lgamma(alpha + count2.at(k).at(0)) - lgamma(alpha);
    }
    score += lgamma(r * alpha) - lgamma(r * alpha + n_i);
  }
  else {
    //have parents
    vector<vector<int>> count;
    count = state_count(data, node_x_vec, parents, n_states);
    int q = count.at(0).size(); //number of states of parents
    int r = count.size(); //number of states of node_x
    vector<float> n_ij(q, 0.0);
    ////cout<<"independent_score_x"<<endl;
    for (int k = 0; k < r; k++) {
      for(int j = 0; j < q; j++) {
        n_ij.at(j) += count.at(k).at(j);
        ////cout<<count.at(k).at(j)<<", ";
      }
      ////cout<<endl;
    }
    for (int j = 0; j < q; j++) { //for each state of parents
      for (int k = 0; k < r; k++) { //for each state of node_x
        score += lgamma(alpha + count.at(k).at(j)) - lgamma(alpha);
        if (isnan(score)) {
          ////////////cout <<ESS / (r * q) + count.at(k).at(j)<< "score is nan" << endl;
        }
      }
      score += lgamma(r * alpha) - lgamma(r * alpha + n_ij.at(j));
    }
    vector<vector<int>> count2;
    count2 = state_count(data, node_y_vec, parents, n_states);
    q = count2.at(0).size(); //number of states of parents
    r = count2.size(); //number of states of node_x
    vector<float> n_ij2(q, 0.0);
    ////cout << "independent_score_y"<<endl;
    for (int k = 0; k < r; k++) {
      for(int j = 0; j < q; j++) {
        n_ij2.at(j) += count2.at(k).at(j);
        ////cout<<count2.at(k).at(j)<<", ";
      }
      ////cout<<endl;
    }
    for (int j = 0; j < q; j++) { //for each state of parents
      for (int k = 0; k < r; k++) { //for each state of node_x
        score += lgamma(alpha + count2.at(k).at(j)) - lgamma(alpha);
        if (isnan(score)) {
          ////////////cout <<ESS / (r * q) + count2.at(k).at(j)<< "score is nan" << endl;
        }
      }
      score += lgamma(r * alpha) - lgamma(r * alpha + n_ij2.at(j));
    }
  }
    //calculate the score
  return score;
}

float natori_dependent_score(const vector<vector<int>> &data, int &node_x, int &node_y, vector<int> &parents, vector<int> &n_states, float &ESS) {
  //return log of the BDeu score of X, Y | Z
  double score = 0.0;
  double alpha = 0.5;

  
  vector<int> children{node_x, node_y};
  if (parents.empty()) {
    //no parents
    vector<vector<int>> count;
    count = state_count(data, children, parents, n_states);
    int r = count.size(); //number of states of node_x]
    if (ESS > -1){
    alpha = ESS / r;
    }
    int n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count.at(k).at(0);
    }
    for (int k = 0; k < r; k++) { //for each state of node_x
        score += lgamma(alpha + count.at(k).at(0)) - lgamma(alpha);
    }
    score += lgamma(r * alpha) - lgamma(r * alpha + n_i);
  }
  else {
    //have parents
    vector<vector<int>> count;
    count = state_count(data, children, parents, n_states);
    int q = count.at(0).size(); //number of states of parents
    int r = count.size(); //number of states of node_x
    if (ESS > -1){
    alpha = ESS / (r * q);
    }
    vector<float> n_ij(q, 0.0);
    ////cout<<"dependent_score"<<endl;
    for (int k = 0; k < r; k++) {
      for(int j = 0; j < q; j++) {
        n_ij.at(j) += count.at(k).at(j);
        ////cout<<count.at(k).at(j)<<", ";
      }
      ////cout<<endl;
    }
    for (int j = 0; j < q; j++) { //for each state of parents
      for (int k = 0; k < r; k++) { //for each state of node_x
        score += lgamma(alpha + count.at(k).at(j)) - lgamma(alpha);
        if (isnan(score)) {
          ////////////cout <<ESS / (r * q) + count.at(k).at(j)<< "score is nan" << endl;
        }
      }
      score += lgamma(r * alpha) - lgamma(r * alpha + n_ij.at(j));
    }
  }
    //calculate the score
    //score = score * 1.000001;
  return score;
}

bool ci_test(const vector<vector<int>> &data, int &node_x, int &node_y, vector<int> &Z, vector<int> &n_states, float &ESS, int &parallel) {
  //CI test for X _|_ Y | Z
  float independent_score = 0.0;
  float dependent_score = 0.0;
  independent_score += natori_independent_score(data, node_x, node_y, Z, n_states, ESS);
  dependent_score += natori_dependent_score(data, node_x, node_y, Z, n_states, ESS);
  if(independent_score > dependent_score){
    cout<< "CI independent:" <<node_x<<" _|_"<<node_y<<" | "<<independent_score<<">"<<dependent_score<< endl;
    return true;
  }else{
    cout<< "CI dependent:" <<node_x<<" _|_"<<node_y<<" | "<<independent_score<<"<"<<dependent_score<< endl;
    return false;
  }
}


double qnorm(double u)
{
	static double a[9] = {	 1.24818987e-4, -1.075204047e-3, 5.198775019e-3,
							-0.019198292004, 0.059054035642,-0.151968751364,
							 0.319152932694,-0.5319230073,   0.797884560593};
	static double b[15] = {	-4.5255659e-5,   1.5252929e-4,  -1.9538132e-5,
							-6.76904986e-4,  1.390604284e-3,-7.9462082e-4,
							-2.034254874e-3, 6.549791214e-3,-0.010557625006,
							 0.011630447319,-9.279453341e-3, 5.353579108e-3,
							-2.141268741e-3, 5.35310549e-4,  0.999936657524};
	double w, y, z;
	int i;

	if(u == 0.)	return 0.5;
	y = u / 2.;
	if(y < -6.)	return 0.;
	if(y > 6.)		return 1.;
	if(y < 0.)		y = - y;
	if(y < 1.)
	{
		w = y * y;
		z = a[0];
		for(i = 1; i < 9; i++)		z = z * w + a[i];
		z *= (y * 2.);
	}
	else
	{
		y -= 2.;
		z = b[0];
		for(i = 1; i < 15; i++)	z = z * y + b[i];
	}

	if(u < 0.)	return (1. - z) / 2.;
	return (1. + z) / 2.;
}

double qchi(double x2, int n)
{
	double w, pw, x, qc;
	int i, i1;

	if(n < 1)
	{
		fprintf(stderr,"Error : 自由度 < 1 in qchi()!\n");
		return 0.;
	}
	if(x2 <= 0.)	return 1.;
	if(x2 > 400.)	return 0.;
	if(n > 10)
	{
		w = 2. / (9. * (double)n);
		pw = pow(x2 / (double)n, 1. / 3.);
		return 1. - qnorm((pw - 1. + w) / sqrt(w));
	}
	w = exp(-x2 / 2.);
	if(n == 2)	return w;
	x = sqrt(x2);
	if(n == 1)	return 2. * (1. - qnorm(x));
	if((n % 2) == 0)
	{
		i1 = 2;
		qc = w;
	}
	else
	{
 		i1 = 1;
		qc = 2. * (1. - qnorm(x));
		w *= (0.797884560750774 / x);
	}
	for(i = i1; i <= n - 2; i += 2)
	{
		w *= (x2 / (double)i);
		qc += w;
	}
	return qc;
}

bool ci_test_for_vstructure_detect_by_Chi_squared_test(const vector<vector<int>> &data, int &node_x, int &node_y, int &node_z, vector<int> &n_states) {
  //CI test for X _|_ Y | Z using Chi-squared test; Test statistics G = 2 * sum_x,y,z(n_xyz * log(n_xyz / (n_xz n_yz / n))) ~ chi^2_{d.f= (|X| - 1) * (|Y| - 1) * |Z|} (from http://www.ai.lab.uec.ac.jp/wp-content/uploads/2019/04/41f06ccbd0ac30d15c7728117770b105.pdf)
  float threthold = 0.05;
  vector<int> children(2);
  children.at(0) = node_x;
  children.at(1) = node_y;
  vector<int> parents(1, node_z);
  vector<vector<int>>count;
  count = state_count(data, children, parents, n_states);
  int r_x = n_states.at(node_x);
  int r_y = n_states.at(node_y);
  int r_z = n_states.at(node_z);
  int n = data.size();
  vector<vector<double>> n_xz(r_x, vector<double>(r_z, 0));
  vector<vector<double>> n_yz(r_y, vector<double>(r_z, 0));
  for (int i = 0; i < r_x; i++) {
    for (int j = 0; j < r_y; j++) {
      for (int k = 0; k < r_z; k++) {
        n_xz.at(i).at(k) += count.at(i * r_y + j).at(k);
        n_yz.at(j).at(k) += count.at(i * r_y + j).at(k);
        //////cout << count.at(i * r_y + j).at(k)<< ", ";
      }
      //////cout << endl;
    }
  }
  vector<double> sum_x(r_z, 0);
  vector<double> sum_y(r_z, 0);
  vector<double> sum_z(r_z, 0);
  for (int i = 0; i < r_x; i++) {
    for (int k = 0; k < r_z; k++) {
      sum_x.at(k) += n_xz.at(i).at(k);
    }
  }
  for (int j = 0; j < r_y; j++) {
    for (int k = 0; k < r_z; k++) {
      sum_y.at(k) += n_yz.at(j).at(k);
    }
  }
  for (int i = 0; i < r_x; i++) {
    for (int j = 0; j < r_y; j++) {
      for (int k = 0; k < r_z; k++) {
        sum_z.at(k) += count.at(i * r_y + j).at(k);
      }
    }
  }
  //////cout<< "n_xz"<<endl;
  for (int i = 0; i < r_x; i++) {
    for (int k = 0; k < r_z; k++) {
      //////cout << n_xz.at(i).at(k)/ sum_x.at(k)<< ", ";
    }
    //////cout << endl;
  }
  //////cout<< "n_yz"<<endl;
  for (int j = 0; j < r_y; j++) {
    for (int k = 0; k < r_z; k++) {
      //////cout << n_yz.at(j).at(k)/ sum_y.at(k)<< ", ";
    }
    //////cout << endl;
  }
  double G = 0.0;
  for (int i = 0; i < r_x; i++) {
    for (int j = 0; j < r_y; j++) {
      for (int k = 0; k < r_z; k++) {
        //////cout << "ijk:" << (double)count.at(i * r_y + j).at(k) /(double)n << "xzyz:" <<((double)n_xz.at(i).at(k) * (double)n_yz.at(j).at(k) / (sum_x.at(k) * sum_y.at(k))) *sum_z.at(k)/(double)n << "ijk/xzyz:"<< ((double)count.at(i * r_y + j).at(k) / (double)n )/ (((double)n_xz.at(i).at(k) * (double)n_yz.at(j).at(k) / (sum_x.at(k) * sum_y.at(k)))*sum_z.at(k)/(double)n)<< endl;
        // G = G + 2 * (double)count.at(i * r_y + j).at(k) * log(((double)count.at(i * r_y + j).at(k) / (double)n )/ (((double)n_xz.at(i).at(k) * (double)n_yz.at(j).at(k) / (sum_x.at(k) * sum_y.at(k)))*sum_z.at(k)/(double)n));
        G = G + 2 * (double)count.at(i * r_y + j).at(k) * log((double)count.at(i * r_y + j).at(k) * (double)n / ((double)n_xz.at(i).at(k) * (double)n_yz.at(j).at(k)));
      }
    }
  }
  //////cout<< "G:" <<G<< endl;
  int dim = (r_x - 1) * (r_y - 1) * r_z;
  double p_value = qchi(G, dim);
  //cout<< "p_value:" <<p_value<< endl;
  bool flag = false;
  if (p_value <= threthold || 1 - p_value <= threthold){
    flag = true;
  }
  if (flag){
    return true;//dependent, reject null hypothesis, find V-structure
    ////cout<< "CI dependent:" <<node_x<<" _|_"<<node_y<<" | "<< node_z << endl;
  }
  else{
    ////cout<< "CI independent:" <<node_x<<" _|_"<<node_y<<" | "<< node_z << endl;
    return false;
  }
}

float localBDeuscore(const vector<vector<int>> &data, int &node_x, vector<int> &parents, vector<int> &n_states) {
  //return log of the BDeu score of X, Y | Z
  double score = 0.0;
  float alpha = 0.5;
  if (parents.empty()) {
    //no parents
    vector<vector<int>> count;
    vector<int> node_x_vec(1, node_x);
    count = state_count(data, node_x_vec, parents, n_states);
    int r = count.size(); //number of states of node_x
    int n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count.at(k).at(0);
    }
    for (int k = 0; k < r; k++) { //for each state of node_x
        score += lgamma(alpha + count.at(k).at(0)) - lgamma(alpha);
        //cout << count.at(k).at(0) << ", ";
        //cout << endl;
    }
    score += lgamma(alpha * r) - lgamma(alpha * r + n_i);
  }
  else {
    //have parents
    vector<vector<int>> count;
    vector<int> node_x_vec(1, node_x);
    count = state_count(data, node_x_vec, parents, n_states);
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
        score += lgamma(alpha + count.at(k).at(j)) - lgamma(alpha);
        if (isnan(score)) {
          ////cout <<ESS / (r * q) + count.at(k).at(j)<< "score is nan" << endl;
        }
        //cout << count.at(k).at(j) << ", ";
      }
      score += lgamma(alpha * r) - lgamma(alpha * r + n_ij.at(j));
      //cout << endl;
    }
  }
    //calculate the score
  return score;
}

bool ci_test_for_vstructure_detect_by_beyesfactor(const vector<vector<int>> &data, int &node_x, int &node_y, int &node_z, vector<int> &n_states) {
  //CI test for X _|_ Y | Z compare X->Z<-Y and X->Z->Y
  float independent_score = 0.0; //X->Z->Y
  float dependent_score = 0.0; //X->Z<-Y

  double alpha = 0.5;

  vector<vector<int>> count_independent;
  vector<int> children_independent = {node_x, node_y, node_z};
  vector<int> parents_independent;
  count_independent = state_count(data, children_independent, parents_independent, n_states);
  int r = count_independent.size(); //number of states of node_x
  int n_i = 0;
  for (int k = 0; k < r; k++) {
    n_i += count_independent.at(k).at(0);
  }
  for (int k = 0; k < r; k++) { //for each state of node_x
      independent_score += lgamma(alpha + count_independent.at(k).at(0)) - lgamma(alpha);
  }
  independent_score += lgamma(alpha) - lgamma(alpha + n_i);
  

  vector<vector<int>> count_dependent_z;
  vector<int> parents_dependent_z = {node_x, node_y};
  vector<int> parents_dependent_x;
  vector<int> parents_dependent_y;
  vector<int> zz = {node_z};
  dependent_score += localBDeuscore(data, node_z, parents_dependent_z, n_states);
  dependent_score += localBDeuscore(data, node_x, parents_dependent_x, n_states);
  //cout<< "xlocal"<<localBDeuscore(data, node_x, parents_dependent_x, n_states)<<endl;
  dependent_score += localBDeuscore(data, node_y, parents_dependent_y, n_states);
  //cout<< "ylocal"<<localBDeuscore(data, node_y, parents_dependent_y, n_states)<<endl;
  // float ess = -2.0;
  // dependent_score += natori_dependent_score(data, node_x, node_y, zz, n_states, ess);
  if(independent_score > dependent_score){
    //cout<< "CI independent:" <<node_x<<" _|_"<<node_y<<" | "<<independent_score<<">"<<dependent_score<< endl;
    return false;
  }else{
    //cout<< "CI dependent:" <<node_x<<" _|_"<<node_y<<" | "<<independent_score<<"<"<<dependent_score<< endl;
    return true;
  }
}


bool ci_test_for_vstructure_detect_by_CMI(const vector<vector<int>> &data, int &node_x, int &node_y, int &node_z, vector<int> &n_states) {
  //CI test for X _|_ Y | Z using CMI; Test statistics cmi(x,y|z) = sum_x,y,z p(x,y,z) * log(p(x,y|z)/p(x|z)p(y|z))
  double threthold = 0.003;
  vector<int> children(2);
  children.at(0) = node_x;
  children.at(1) = node_y;
  vector<int> parents(1, node_z);
  vector<vector<int>>count;
  count = state_count(data, children, parents, n_states);
  int r_x = n_states.at(node_x);
  int r_y = n_states.at(node_y);
  int r_z = n_states.at(node_z);
  int n = data.size();
  vector<vector<double>> n_xz(r_x, vector<double>(r_z, 0));
  vector<vector<double>> n_yz(r_y, vector<double>(r_z, 0));
  for (int i = 0; i < r_x; i++) {
    for (int j = 0; j < r_y; j++) {
      for (int k = 0; k < r_z; k++) {
        n_xz.at(i).at(k) += count.at(i * r_y + j).at(k);
        n_yz.at(j).at(k) += count.at(i * r_y + j).at(k);
        //////cout << count.at(i * r_y + j).at(k)<< ", ";
      }
      //////cout << endl;
    }
  }
  vector<double> sum_x(r_z, 0);
  vector<double> sum_y(r_z, 0);
  vector<double> sum_z(r_z, 0);
  for (int i = 0; i < r_x; i++) {
    for (int k = 0; k < r_z; k++) {
      sum_x.at(k) += n_xz.at(i).at(k);
    }
  }
  for (int j = 0; j < r_y; j++) {
    for (int k = 0; k < r_z; k++) {
      sum_y.at(k) += n_yz.at(j).at(k);
    }
  }
  for (int i = 0; i < r_x; i++) {
    for (int j = 0; j < r_y; j++) {
      for (int k = 0; k < r_z; k++) {
        sum_z.at(k) += count.at(i * r_y + j).at(k);
      }
    }
  }


  vector<vector<double>> p_xyz(r_x * r_y, vector<double>(r_z, 0.0));
  vector<vector<double>> p_xy_z(r_x * r_y, vector<double>(r_z, 0.0));
  vector<vector<double>> p_xz(r_x, vector<double>(r_z, 0.0));
  vector<vector<double>> p_yz(r_y, vector<double>(r_z, 0.0));
  //cout << "p_xy_z" << endl;
  for (int i = 0; i < r_x; i++) {
    for (int j = 0; j < r_y; j++) {
      for (int k = 0; k < r_z; k++) {
        p_xyz.at(i * r_y + j).at(k) = (double)count.at(i * r_y + j).at(k) / (double)n; 
        p_xy_z.at(i * r_y + j).at(k) = (double)count.at(i * r_y + j).at(k) / (double)sum_z.at(k);
        //cout << p_xy_z.at(i * r_y + j).at(k) << ", "; 
      }
      //cout << endl;
    }
  }
  //cout << "p_xz" << endl;
  for (int i = 0; i < r_x; i++) {
    for (int k = 0; k < r_z; k++) {
      p_xz.at(i).at(k) = (double)n_xz.at(i).at(k) / (double)sum_x.at(k);
      //cout << p_xz.at(i).at(k) << ", ";
    }
    //cout << endl;
  }
  //cout << "p_yz" << endl;
  for (int j = 0; j < r_y; j++) {
    for (int k = 0; k < r_z; k++) {
      p_yz.at(j).at(k) = (double)n_yz.at(j).at(k) / (double)sum_y.at(k);
      //cout << p_yz.at(j).at(k) << ", ";
    }
    //cout << endl;
  }
  //cout << "log" << endl;
  double cmi = 0.0;
  for (int i = 0; i < r_x; i++) {
    for (int j = 0; j < r_y; j++) {
      for (int k = 0; k < r_z; k++) {
        double temp = 0;
        temp = p_xyz.at(i * r_y + j).at(k) * log(p_xy_z.at(i * r_y + j).at(k) / (p_xz.at(i).at(k) * p_yz.at(j).at(k)));
        cmi = cmi + temp;
        ////cout <<"(" << p_xy_z.at(i * r_y + j).at(k) / (p_xz.at(i).at(k) * p_yz.at(j).at(k)) <<"," << log(p_xy_z.at(i * r_y + j).at(k) / (p_xz.at(i).at(k) * p_yz.at(j).at(k)))<<")";
        ////cout << temp << ", ";
      }
      ////cout << endl;
    }
  }


  
  //cout << "CMI: " << cmi << endl;
  if(cmi < threthold){
    //cout<< "CI independent:" <<node_x<<" _|_"<<node_y<<" | "<<node_z<< endl;
    return false;
  }else{
    //cout<< "CI dependent:" <<node_x<<" _|_"<<node_y<<" | "<<node_z<< endl;
    return true;
  }
}

bool ci_test_for_vstructure_detect_by_MI(const vector<vector<int>> &data, int &node_x, int &node_y, vector<int> &n_states) {
  //CI test for X _|_ Y using MI; Test statistics cmi(x,y) = sum_x,y,z p(x,y) * log(p(x,y)/p(x)p(y))
  double threthold = 0.003;
  vector<int> children(2);
  children.at(0) = node_x;
  children.at(1) = node_y;
  vector<int> parents;
  vector<vector<int>>count;
  count = state_count(data, children, parents, n_states);
  int r_x = n_states.at(node_x);
  int r_y = n_states.at(node_y);
  int n = data.size();

  vector<double> p_xy(r_x * r_y, 0.0);
  vector<double> p_x(r_x, 0.0);
  vector<double> p_y(r_y, 0.0);
  for (int i = 0; i < r_x; i++) {
    for (int j = 0; j < r_y; j++) {
      p_xy.at(i * r_y + j) = (double)count.at(i * r_y + j).at(0) / (double)n; 
    }
  }
  cout << "p_xy" << endl;
  for (int i = 0; i < r_x; i++) {
    for (int j = 0; j < r_y; j++) {
      p_x.at(i) += p_xy.at(i * r_y + j);
      p_y.at(j) += p_xy.at(i * r_y + j);
      cout << p_xy.at(i * r_y + j) << ", ";
    }
    cout << endl;
  }
  double cmi = 0.0;
  for (int i = 0; i < r_x; i++) {
    for (int j = 0; j < r_y; j++) {
      double temp = 0;
      temp = p_xy.at(i * r_y + j) * log(p_xy.at(i * r_y + j) / (p_x.at(i) * p_y.at(j)));
      cmi = cmi + temp;
    }
  }

  if(cmi < threthold){
    cout<< "CI independent:" <<node_x<<" _|_"<<node_y << cmi<< endl;
    return true;
  }else{
    cout<< "CI dependent:" <<node_x<<" _|_"<<node_y<< cmi<< endl;
    return false;
  }
}

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
            ////////////cout<< "A2_removed:" <<Y<<"->"<<Z<< endl;
          }
        }
      }
    }
  }
  return;
}

void orientation_B2(PDAG &Gall, vector<int> &Gs, vector<vector<bool>> &deletededges, const vector<vector<int>> &data, vector<int> &n_states, float &ESS, int &parallel) {
  /*
      orient edges in a PDAG to a maximally oriented graph.
      orient rules are based on rule 1~3 from Meek,C.:Causal Inference and Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in Artificial Inteligence (UAl-95),p.403-410 (195)
  */
  //for each X-Z-Y (X and Y is not adjecent), find V-structure and orient as X -> Z <- Y
  for (auto& X : Gs) {
    for (auto& Z : Gall.undirected_neighbors(X)) {
      for (auto& Y : Gall.undirected_neighbors(Z)) {
        //if (X != Y && !Gall.has_edge(X, Y) && !Gall.has_edge(Y, X) && Gall.has_edge(X, Z) && Gall.has_edge(Z, X) && Gall.has_edge(Y, Z) && Gall.has_edge(Z, Y)) {
        if (X != Y && !Gall.has_edge(X, Y) && !Gall.has_edge(Y, X) && Gall.has_edge(X, Z) && Gall.has_edge(Y, Z) && (Gall.has_edge(Z, Y) || Gall.has_edge(Z, X))) {
          vector<int> z = {Z};
          vector<int> v;
          cout<< "V-structure think:" <<X<<"->"<<Z<<"<-"<<Y<< endl;
          //if(!ci_test(data, X, Y, z, n_states, ESS)){
          if(ci_test(data, X, Y, v, n_states, ESS, parallel)){
          //if(ci_test_for_vstructure_detect_by_MI(data, X, Y, n_states)){
          //if(ci_test_for_vstructure_detect_by_CMI(data, X, Y, v, n_states)){
          //if(ci_test_for_vstructure_detect_by_Chi_squared_test(data, X, Y, Z, n_states)){
          //if(ci_test_for_vstructure_detect_by_CMI(data, X, Y, Z, n_states)){ 
          //if(ci_test_for_vstructure_detect_by_beyesfactor(data, X, Y, Z, n_states)){
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
        if (!Gall.has_edge(Y, X) && Gall.has_edge(X, Y)) {
          for (auto& Z : Gall.undirected_neighbors(Y)) {
            if (!Gall.has_edge(X, Z) && !Gall.has_edge(Z, X) && Z != X) {
              Gall.remove_edge(Z, Y);
              cout<< "R1:" <<Y<<"->"<<Z<<"|"<<X<< endl;
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
                cout<< "R3:" <<Y<<"->"<<W<<"|"<<"X="<<X<<"Z="<<Z<< endl;
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

void recursive_search(const vector<vector<int>> &data, PDAG &Gall, vector<int> Gs, vector<int> &Gex, int N, vector<int> &n_states, float &ESS, int &parallel) {
  int n_node = data.at(0).size();
  //print
  ////cout<< "N=" <<N<< ", ";
  ////cout<< "Gs: ";
  for (auto& node : Gs) {
    ////cout<< node << ", ";
  }
  ////cout<< "Gex: ";
  for (auto& node : Gex) {
    ////cout<< node << ", ";
  }
  ////cout<< endl;

  //stage 0: exit condition
  bool exitcondition = true;
  for (int i = 0; i < Gs.size(); i++) {
    if (Gall.predecessors(Gs.at(i)).size() > N){
      exitcondition = false;
      break;
    }
  }
  if (exitcondition) {
    return;
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
                  if (!deletededges.at(node_x).at(node_y)){
                    if (ci_test(data, node_x, node_y, selected_z, n_states, ESS, parallel)) {
                      Gall.remove_edge(node_x, node_y);
                      deletededges.at(node_x).at(node_y) = true;
                      deletededges.at(node_y).at(node_x) = true;
                      cout<< "A1_removed:" <<node_x<<"-"<<node_y<< endl;
                      //transive_cut();
                    }
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
    for (auto& node_x : Gs) {
      if (!deletededges.at(node_x).at(node_y) && node_x != node_y && (Gall.has_edge(node_x, node_y) || Gall.has_edge(node_y, node_x))) {
        if (N == 0) {
          vector<int> S;
            if (ci_test(data, node_x, node_y, S, n_states, ESS, parallel)) {
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
              if (ci_test(data, node_x, node_y, selected_z, n_states, ESS, parallel)) {
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
  }
  
  //stage B2: orient edges in Gs using orientation rules R1~R3
  orientation_B2(Gall, Gs, deletededges, data, n_states, ESS, parallel);
  //stage B3: Group the nodes having the lowest topological order into a descendant substructure Gd
  vector<int> Gd;
  vector<vector<int>> g_ex_connection;
  order_grouping(Gall, Gs, Gd, g_ex_connection);
  //stage C: Ancestor sub-structure decomposition
  vector<int> Gexd;
  for (int i = 0; i < g_ex_connection.size(); i++) {
    for (int j = 0; j< g_ex_connection.at(i).size(); j++) {
      Gexd.push_back(g_ex_connection.at(i).at(j));
    }
    recursive_search(data, Gall, g_ex_connection.at(i), Gex, N + 1, n_states, ESS);
  }
  sort(Gexd.begin(), Gexd.end());
  for (int i = 0; i < Gexd.size(); i++) {
    bool check = true;
    for(int j = 0; j < Gex.size(); j++){
      if (Gexd.at(i) == Gex.at(j)){
        check = false;
        break;
      }
    }
    if (check){
      Gex.push_back(Gexd.at(i));
    } 
  }


  //stage D:  Descendantsub-structure decomposition
  recursive_search(data, Gall, Gd, Gex, N + 1, n_states, ESS);
  return;
}

py::array_t<bool> RAI(py::array_t<int> data, py::array_t<int> n_states, float ESS, int parallel) {
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
  recursive_search(data_vec, Gall, Gs, Gex, 0, n_states_vec, ESS, parallel);

  //translate Gall to py::array_t (this is not optimal but I don't know how to use pybind11::array_t)
  for (size_t i = 0; i < n_node; i++) {
    for (size_t j = 0; j < n_node; j++) {
      prt_endg[i * n_node + j] = Gall.g.at(i).at(j);
    }
  }
  return endg;
}

