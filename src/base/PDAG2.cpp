#include "base/PDAG2.h"

#include <iostream>
#include <queue>

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

set<int> PDAG::successors(int i) { return successor_sets[i]; }

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

void PDAG::remove_edge(int i, int j) {
  if (!g.at(i).at(j)) return;
  g.at(i).at(j) = false;
  if (g.at(j).at(i)) {
    successor_sets[j].insert(i);
  } else {
    successor_sets[i].erase(j);
  }
}

void PDAG::remove_edge_completedly(int i, int j) {
  g.at(i).at(j) = false;
  g.at(j).at(i) = false;
}

void PDAG::add_edge(int i, int j) {
  if (g.at(i).at(j)) return;
  g.at(i).at(j) = true;
  if (g.at(j).at(i)) {
    successor_sets[j].erase(i);
  } else {
    successor_sets[i].insert(j);
  }
}

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

bool PDAG::has_cycle() {
  int n = g.size();
  vector<int> indeg(n);
  for (int i = 0; i < n; i++) {
    for (int j : successors(i)) {
      indeg[j]++;
    }
  }
  queue<int> que;
  for (int i = 0; i < n; i++) {
    if (!indeg[i]) que.push(i);
  }
  int cnt = 0;
  while (!que.empty()) {
    int v = que.front();
    que.pop();
    cnt++;
    for (int u : successors(v)) {
      indeg[u]--;
      if (!indeg[u]) que.push(u);
    }
  }
  return cnt != n;
}