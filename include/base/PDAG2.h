#pragma once
#include <vector>
using namespace std;

struct PDAG {
  vector<vector<bool>> g;
  // コンストラクタ
  PDAG();
  // コピーコンストラクタ
  PDAG(const PDAG &old);
  // 代入演算子
  PDAG &operator=(const PDAG &a);
  // デストラクタ
  ~PDAG() = default;

  // return the list of successors of node i (include undirected edge)
  vector<int> successors(int i);

  // return the list of predecessors of node i (include undirected edge)
  vector<int> predecessors(int i);

  // return the list of neighbors {j} of node i (j -> i or i -> j)
  vector<int> neighbors(int i);

  // return the list of undirected neighbors of node i
  vector<int> undirected_neighbors(int i);

  // remove the edge i -> j
  void remove_edge(int i, int j);

  // remove edge between i and j
  void remove_edge_completedly(int i, int j);

  void add_edge(int i, int j);

  bool has_edge(int i, int j);
  bool has_directed_edge(int i, int j);
  bool has_undirected_edge(int i, int j);

  // check if there is a directed path from X to Y using DFS
  bool has_directed_path(int X, int Y);
};