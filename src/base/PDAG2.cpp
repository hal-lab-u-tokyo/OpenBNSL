#include "base/PDAG2.h"

PDAG::PDAG() {
  // cout << "normal constructor called" << endl;
}

PDAG::PDAG(const PDAG& old) {
  // cout << "copy constructor called" << endl;
  g = old.g;
}

PDAG& PDAG::operator=(const PDAG& a) {
  if (this != &a) g = a.g;
  return *this;
}

vector<int> PDAG::successors(int i) {
  vector<int> succ;
  for (int j = 0; j < (int)g.size(); j++) {
    if (g.at(i).at(j)) {
      succ.push_back(j);
    }
  }
  return succ;
}

vector<int> PDAG::predecessors(int i) {
  vector<int> pred;
  for (int j = 0; j < (int)g.size(); j++) {
    if (g.at(j).at(i)) {
      pred.push_back(j);
    }
  }
  return pred;
}

vector<int> PDAG::neighbors(int i) {
  vector<int> neigh;
  for (int j = 0; j < (int)g.size(); j++) {
    if (g.at(j).at(i) || g.at(i).at(j)) {
      neigh.push_back(j);
    }
  }
  return neigh;
}

vector<int> PDAG::undirected_neighbors(int i) {
  vector<int> neigh;
  for (int j = 0; j < (int)g.size(); j++) {
    if (g.at(i).at(j) && g.at(j).at(i)) {
      neigh.push_back(j);
    }
  }
  return neigh;
}

void PDAG::remove_edge(int i, int j) { g.at(i).at(j) = false; }

void PDAG::remove_edge_completedly(int i, int j) {
  g.at(i).at(j) = false;
  g.at(j).at(i) = false;
}

void PDAG::add_edge(int i, int j) { g.at(i).at(j) = true; }

bool PDAG::has_edge(int i, int j) { return g.at(i).at(j); }

bool PDAG::has_directed_edge(int i, int j) {
  if (g.at(i).at(j) && !g.at(j).at(i)) {
    return true;
  }
  return false;
}

bool PDAG::has_undirected_edge(int i, int j) {
  return g.at(i).at(j) && g.at(j).at(i);
}

bool PDAG::has_directed_path(int X, int Y) {
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
    for (auto& succ : successors(node)) {
      if (visited.at(succ) == 0 && has_directed_edge(node, succ)) {
        stack.push_back(succ);
      }
    }
  }
  return false;
}