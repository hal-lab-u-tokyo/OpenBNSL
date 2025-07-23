#pragma once
using namespace std;

#include <vector>

struct PDAG_nishikori {
  vector<vector<bool>> g;
  // コンストラクタ
  PDAG_nishikori() {}
  // コピーコンストラクタ
  PDAG_nishikori(const PDAG_nishikori &old) { g = old.g; }
  // 代入演算子
  PDAG_nishikori &operator=(const PDAG_nishikori &a) {
    if (this != &a) g = a.g;
    return *this;
  }
  // デストラクタ
  ~PDAG_nishikori() = default;

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

  bool has_path(
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
        if (visited.at(succ) == 0 && has_edge(node, succ)) {
          stack.push_back(succ);
        }
      }
    }
    return false;
  }

  bool has_connection(
      int X,
      int Y) {  // check if there is connection between X and Y using DFS
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
        if (visited.at(succ) == 0 &&
            (has_edge(node, succ) || has_edge(succ, node))) {
          stack.push_back(succ);
        }
      }
    }
    return false;
  }
};