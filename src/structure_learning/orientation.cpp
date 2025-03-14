#include "structure_learning/orientation.h"

#include "structure_learning/constants.h"

void orientation(PDAG &G, const vector<int> &sepsets) {
  // for each X-Z-Y (X and Y is not adjecent), find V-structure and orient as
  // X
  // -> Z <- Y
  int n_node = G.g.size();
  auto v_structure = vector<vector<bool>>(n_node, vector<bool>(n_node));
  for (int X = 0; X < n_node; X++) {
    for (int Z : G.undirected_neighbors(X)) {
      for (int Y : G.undirected_neighbors(Z)) {
        if (X == Y || G.has_edge(X, Y) || G.has_edge(Y, X)) continue;
        bool in_sepset = false;
        for (int i = 0; i < max_level; i++) {
          int XYmin = (X < Y ? X : Y);
          int XYmax = (X < Y ? Y : X);
          int id = sepsets[(XYmin * n_node + XYmax) * max_level + i];
          if (id == -1) break;
          if (id == Z) {
            in_sepset = true;
            break;
          }
        }
        if (!in_sepset) {
          v_structure[X][Z] = true;
          v_structure[Y][Z] = true;
          // cout << "V-structure found:" << X << "->" << Z << "<-" << Y <<
          // endl;
        }
      }
    }
  }
  for (int i = 0; i < n_node; i++) {
    for (int j = 0; j < n_node; j++) {
      if (v_structure[i][j] && !v_structure[j][i]) {
        G.remove_edge(j, i);
      }
    }
  }
  bool flag = true;
  while (flag) {
    flag = false;
    // Rule 1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
    for (int X = 0; X < n_node; X++) {
      for (int Y : G.successors(X)) {
        if (!G.has_directed_edge(X, Y)) continue;
        for (int Z : G.undirected_neighbors(Y)) {
          if (!G.has_edge(X, Z) && !G.has_edge(Z, X) && Z != X) {
            G.remove_edge(Z, Y);
            // cout << "R1:" << Y << "->" << Z << endl;
            flag = true;
          }
        }
      }
    }
    // Rule 2: X - Y and if there is a directed path from X to Y, then X -> Y
    for (int X = 0; X < n_node; X++) {
      for (int Y : G.undirected_neighbors(X)) {
        if (G.has_directed_path(X, Y)) {
          G.remove_edge(Y, X);
          // cout << "R2:" << X << "->" << Y << endl;
          flag = true;
        }
      }
    }
    // Rule 3: for each X->W<-Z X-Y-Z Y-W, orient Y->W
    for (int X = 0; X < n_node; X++) {
      for (int Y : G.undirected_neighbors(X)) {
        for (int Z : G.undirected_neighbors(Y)) {
          if (Z == X || G.has_edge(X, Z) || G.has_edge(Z, X)) continue;
          // X-Y-Z
          for (int W : G.undirected_neighbors(Y)) {
            if (W != X && W != Z && G.has_directed_edge(X, W) &&
                G.has_directed_edge(Z, W)) {
              G.remove_edge(W, Y);
              // cout << "R3:" << Y << "->" << W << endl;
              flag = true;
            }
          }
        }
      }
    }
  }
  return;
}

void orientation(int level, PDAG &G, const vector<int> &sepsets) {
  // for each X-Z-Y (X and Y is not adjecent), find V-structure and orient as
  // X
  // -> Z <- Y
  int n_node = G.g.size();
  auto v_structure = vector<vector<bool>>(n_node, vector<bool>(n_node));
  for (int X = 0; X < n_node; X++) {
    for (int Y = X + 1; Y < n_node; Y++) {
      if (G.has_edge(X, Y) || G.has_edge(Y, X)) continue;
      int sepset_level = 0;
      for (int i = 0; i < max_level; i++) {
        int id = sepsets[(X * n_node + Y) * max_level + i];
        if (id == -1) break;
        sepset_level++;
      }
      if (sepset_level != level) continue;
      for (int Z : G.undirected_neighbors(X)) {
        if (!G.has_undirected_edge(Y, Z)) continue;
        bool in_sepset = false;
        for (int i = 0; i < level; i++) {
          int id = sepsets[(X * n_node + Y) * max_level + i];
          if (id == Z) {
            in_sepset = true;
            break;
          }
        }
        if (!in_sepset) {
          v_structure[X][Z] = true;
          v_structure[Y][Z] = true;
          // cout << "V-structure found:" << X << "->" << Z << "<-" << Y <<
          // endl;
        }
      }
    }
  }
  for (int i = 0; i < n_node; i++) {
    for (int j = 0; j < n_node; j++) {
      if (v_structure[i][j] && !v_structure[j][i]) {
        G.remove_edge(j, i);
      }
    }
  }
  bool flag = true;
  while (flag) {
    flag = false;
    // Rule 1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
    for (int X = 0; X < n_node; X++) {
      for (int Y : G.successors(X)) {
        if (!G.has_directed_edge(X, Y)) continue;
        for (int Z : G.undirected_neighbors(Y)) {
          if (!G.has_edge(X, Z) && !G.has_edge(Z, X) && Z != X) {
            G.remove_edge(Z, Y);
            // cout << "R1:" << Y << "->" << Z << endl;
            flag = true;
          }
        }
      }
    }
    // Rule 2: X - Y and if there is a directed path from X to Y, then X -> Y
    for (int X = 0; X < n_node; X++) {
      for (int Y : G.undirected_neighbors(X)) {
        if (G.has_directed_path(X, Y)) {
          G.remove_edge(Y, X);
          // cout << "R2:" << X << "->" << Y << endl;
          flag = true;
        }
      }
    }
    // Rule 3: for each X->W<-Z X-Y-Z Y-W, orient Y->W
    for (int X = 0; X < n_node; X++) {
      for (int Y : G.undirected_neighbors(X)) {
        for (int Z : G.undirected_neighbors(Y)) {
          if (Z == X || G.has_edge(X, Z) || G.has_edge(Z, X)) continue;
          // X-Y-Z
          for (int W : G.undirected_neighbors(Y)) {
            if (W != X && W != Z && G.has_directed_edge(X, W) &&
                G.has_directed_edge(Z, W)) {
              G.remove_edge(W, Y);
              // cout << "R3:" << Y << "->" << W << endl;
              flag = true;
            }
          }
        }
      }
    }
  }
  return;
}
