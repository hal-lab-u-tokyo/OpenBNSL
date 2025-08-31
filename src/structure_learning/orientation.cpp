#include "structure_learning/orientation.h"

#include <cassert>
#include <cstddef>
#include <iostream>

#include "structure_learning/constants.h"

void orientation(PDAG &G, const vector<int> &pairs,
                 const vector<int> &sepsets) {
  // for each X-Z-Y (X and Y is not adjecent), find V-structure and orient as
  // X
  // -> Z <- Y
  int n_node = G.g.size();
  auto v_structure = vector<vector<bool>>(n_node, vector<bool>(n_node));
  for (size_t pair_idx = 0; pair_idx < pairs.size() / 2; pair_idx++) {
    int X = pairs[pair_idx * 2], Y = pairs[pair_idx * 2 + 1];
    for (int Z : G.neighbors(X)) {
      if (!G.has_edge(X, Z) || !G.has_edge(Y, Z)) {
        continue;
      }
      int sep_all = sepsets[pair_idx * (n_node + 1)];
      int sep_cnt = sepsets[pair_idx * (n_node + 1) + Z + 1];
      // cout << "X, Z, Y, sep_all, sep_cnt: " << X << ' ' << Z << ' ' << Y << '
      // '
      //      << sep_all << ' ' << sep_cnt << endl;
      if (sep_cnt * 2 <= sep_all) {
        G.remove_edge(Z, X);
        G.remove_edge(Z, Y);
        if (G.has_cycle()) {
          G.add_edge(Z, X);
          G.add_edge(Z, Y);
          // cout << "cycle detected" << endl;
          assert(!G.has_cycle());
        }
      }
    }
    // for (int X = 0; X < n_node; X++) {
    //   for (int Z : G.undirected_neighbors(X)) {
    //     for (int Y : G.undirected_neighbors(Z)) {
    //       if (X == Y || G.has_edge(X, Y) || G.has_edge(Y, X)) continue;
    //       bool in_sepset = false;
    //       for (int i = 0; i < max_level; i++) {
    //         int XYmin = (X < Y ? X : Y);
    //         int XYmax = (X < Y ? Y : X);
    //         int id = sepsets[(XYmin * n_node + XYmax) * max_level + i];
    //         if (id == -1) break;
    //         if (id == Z) {
    //           in_sepset = true;
    //           break;
    //         }
    //       }
    //       if (!in_sepset) {
    //         if (G.has_directed_edge(Z, X) || G.has_directed_edge(Z, Y)) {
    //           cout << "conflict" << endl;
    //         } else {
    //         }
    //         // v_structure[X][Z] = true;
    //         // v_structure[Y][Z] = true;
    //         // cout << "V-structure found:" << X << "->" << Z << "<-" << Y <<
    //         // endl;
    //       }
    //     }
    //   }
    // }
    // for (int i = 0; i < n_node; i++) {
    //   for (int j = 0; j < n_node; j++) {
    //     if (v_structure[i][j] && !v_structure[j][i]) {
    //       G.remove_edge(j, i);
    //       if (G.has_cycle()) {
    //         G.add_edge(j, i);
    //         cout << "cycle detected" << endl;
    //         assert(!G.has_cycle());
    //       }
    //     }
    //   }
    // }
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
            if (G.has_cycle()) {
              G.add_edge(Z, Y);
              // cout << "cycle detected" << endl;
              assert(!G.has_cycle());
            } else {
              // cout << "R1:" << Y << "->" << Z << endl;
              flag = true;
            }
          }
        }
      }
    }
    // Rule 2: X - Y and if there is a directed path from X to Y, then X -> Y
    for (int X = 0; X < n_node; X++) {
      for (int Y : G.undirected_neighbors(X)) {
        if (G.has_directed_path(X, Y)) {
          G.remove_edge(Y, X);
          if (G.has_cycle()) {
            G.add_edge(Y, X);
            // cout << "cycle detected" << endl;
            assert(!G.has_cycle());
          } else {
            // cout << "R2:" << X << "->" << Y << endl;
            flag = true;
          }
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
              if (G.has_cycle()) {
                G.add_edge(W, Y);
                // cout << "cycle detected" << endl;
                assert(!G.has_cycle());
              } else {
                // cout << "R3:" << Y << "->" << W << endl;
                flag = true;
              }
            }
          }
        }
      }
    }
  }
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
