using namespace std;
#include "structure_learning/RAI_nishikori.h"

#include <cassert>
#include <set>
#include <stack>
#include <string>
#include <vector>

#include "base/PDAG_nishikori.h"
// the following includes are for permutation and combination algorithms
#include <algorithm>
#include <functional>
#include <iostream>

// for gamma function
#include <cmath>

// input: data: np.ndarray,  shape: (n: number of variables, d: number of
// samples) output: leard PDAG_nishikori
// Gall.at(i).at(j)==1 means there is an edge i -> j

void dfs(int now, vector<vector<int>> &G, vector<bool> &visited,
         vector<int> &order) {
  visited.at(now) = true;
  for (int i = 0; i < (int)G.at(now).size(); i++) {
    if (!visited.at(G.at(now).at(i))) {
      dfs(G.at(now).at(i), G, visited, order);
    }
  }
  order.push_back(now);
}

void rdfs(int now, vector<vector<int>> &rG, vector<bool> &visited, int k,
          vector<int> &comp) {
  visited.at(now) = true;
  comp.at(now) = k;
  for (int i = 0; i < (int)rG.at(now).size(); i++) {
    if (!visited.at(rG.at(now).at(i))) {
      rdfs(rG.at(now).at(i), rG, visited, k, comp);
    }
  }
}

vector<int> SCC(vector<vector<int>> &G, vector<vector<int>> &rG) {
  int n_node = G.size();
  vector<bool> visited(n_node, false);
  vector<int> comp(n_node, 0);
  vector<int> order;
  order.clear();
  // dfs
  for (int i = 0; i < n_node; i++) {
    if (!visited.at(i)) {
      dfs(i, G, visited, order);
    }
  }
  fill(visited.begin(), visited.end(), false);
  // dfs2
  int k = 0;
  for (int i = n_node - 1; i >= 0; i--) {
    if (!visited.at(order.at(i))) {
      rdfs(order.at(i), rG, visited, k, comp);
      k++;
    }
  }
  return comp;
}  // lowest topological order -> highest number in comp

vector<vector<int>> adjmat2listmat(
    vector<vector<bool>> &adjmat) {  //隣接行列表現to隣接リスト表現
  int n_node = adjmat.size();
  vector<vector<int>> listmat(n_node);
  for (int i = 0; i < n_node; i++) {
    for (int j = 0; j < n_node; j++) {
      if (adjmat.at(i).at(j)) {
        listmat.at(i).push_back(j);  // i to j
      }
    }
  }
  return listmat;
}

vector<vector<int>> adjmat2listmat_reverse(
    vector<vector<bool>> &adjmat) {  //隣接行列表現to隣接リスト表現(reversed)
  int n_node = adjmat.size();
  vector<vector<int>> listmat(n_node);
  for (int i = 0; i < n_node; i++) {
    for (int j = 0; j < n_node; j++) {
      if (adjmat.at(i).at(j)) {
        listmat.at(j).push_back(i);  // reverse edge j to i
      }
    }
  }
  return listmat;
}

vector<vector<bool>> make_gs_graph(vector<vector<bool>> &Gallgraph,
                                   vector<int> &gs) {
  vector<vector<bool>> gsgraph(gs.size(), vector<bool>(gs.size(), false));
  for (int i = 0; i < (int)gs.size(); i++) {
    for (int j = 0; j < (int)gs.size(); j++) {
      gsgraph.at(i).at(j) = Gallgraph.at(gs.at(i)).at(gs.at(j));
    }
  }
  return gsgraph;
}

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

void order_grouping(PDAG_nishikori &Gall, vector<int> &Gs, vector<int> &Gd,
                    vector<vector<int>> &g_ex_connection) {
  vector<vector<bool>> adjmat = make_gs_graph(Gall.g, Gs);
  vector<vector<int>> listmat = adjmat2listmat(adjmat);
  vector<vector<int>> rlistmat = adjmat2listmat_reverse(adjmat);
  vector<int> comp = SCC(listmat, rlistmat);
  PDAG_nishikori Gss;
  Gss.g = make_gs_graph(Gall.g, Gs);
  bool flag = true;
  int order = 0;
  vector<vector<int>> g_subs;
  while (flag) {
    flag = false;
    for (int i = 0; i < (int)comp.size(); i++) {
      if (comp.at(i) == order) {  // node Gs.at(i) is in the order-th group
        if (flag == false) {
          g_subs.push_back(vector<int>());
        }
        g_subs.at(order).push_back(i);  // g_subsはGsのindexを格納
        flag = true;
      }
    }
    order = order + 1;
  }
  for (int i = 0; i < (int)comp.size(); i++) {
  }

  //子孫部分集合の分離
  //各クラスタについてその全てのノードが他のクラスタに子ノードを持たないもの(クラスタ内の一つのノードを持ってきたときに他のクラスタのどれか一つのノードへのpathがないもの)を判定
  vector<bool> is_Gd(g_subs.size(), true);
  for (int i = 0; i < (int)g_subs.size(); i++) {
    for (int j = 0; j < (int)g_subs.size(); j++) {
      if (i != j) {
        if (Gss.has_path(g_subs.at(i).at(0),
                         g_subs.at(j).at(0))) {  // index gss と gall の不一致
          is_Gd.at(i) = false;
          break;
        }
      }
    }
  }
  //それら全てをG_dとする
  vector<int> Gex;
  for (int i = 0; i < (int)g_subs.size(); i++) {
    if (is_Gd.at(i)) {
      for (int j = 0; j < (int)g_subs.at(i).size(); j++) {
        Gd.push_back(g_subs.at(i).at(j));
      }
    } else {
      for (int j = 0; j < (int)g_subs.at(i).size(); j++) {
        Gex.push_back(g_subs.at(i).at(j));
      }
    }
  }
  // G_dが全てなら終了
  bool temp = true;
  for (int i = 0; i < (int)is_Gd.size(); i++) {
    if (!is_Gd.at(i)) {
      temp = false;
      break;
    }
  }
  if (temp) {
    return;
  }

  // G_d以外について，クラスタ0から連結しているクラスタをまとめていく
  PDAG_nishikori Gss_ex;
  Gss_ex.g = make_gs_graph(Gss.g, Gex);
  for (int i = 0; i < (int)Gss_ex.g.size(); i++) {
    for (int j = 0; j < (int)Gss_ex.g.size(); j++) {
      if (Gss_ex.g.at(i).at(j)) {
        Gss_ex.g.at(j).at(i) = true;
      }
    }
  }  // adjmat_ex is undirected graph
  vector<int> Gex_reverseindexlist(Gss.g.size(),
                                   -1);  //.at(Gs_node_index) = Gex_node_index
  for (int i = 0; i < (int)Gex.size(); i++) {
    Gex_reverseindexlist.at(Gex.at(i)) = i;
  }

  vector<bool> flag_ex_vec(g_subs.size(), true);
  for (int i = 0; i < (int)g_subs.size(); i++) {
    if (flag_ex_vec.at(i) && !is_Gd.at(i)) {
      flag_ex_vec.at(i) = false;
      g_ex_connection.push_back(vector<int>());
      for (int l = 0; l < (int)g_subs.at(i).size(); l++) {
        g_ex_connection.back().push_back(g_subs.at(i).at(l));
      }
      for (int j = 0; j < (int)g_subs.size(); j++) {
        if (i != j && !is_Gd.at(j) && !g_subs.at(i).empty() &&
            !g_subs.at(j).empty()) {
          if (Gss_ex.has_connection(
                  Gex_reverseindexlist.at(g_subs.at(i).at(0)),
                  Gex_reverseindexlist.at(g_subs.at(j).at(0)))) {
            for (int k = 0; k < (int)g_subs.at(j).size(); k++) {
              g_ex_connection.back().push_back(g_subs.at(j).at(k));
            }
            flag_ex_vec.at(j) = false;
          }
        }
      }
      sort(g_ex_connection.back().begin(), g_ex_connection.back().end());
    }
  }
  // Gs, g_ex_connectionのindexをGallのindexに変換
  for (int i = 0; i < (int)Gd.size(); i++) {
    Gd.at(i) = Gs.at(Gd.at(i));
  }
  for (int i = 0; i < (int)g_ex_connection.size(); i++) {
    for (int j = 0; j < (int)g_ex_connection.at(i).size(); j++) {
      g_ex_connection.at(i).at(j) = Gs.at(g_ex_connection.at(i).at(j));
    }
  }
  return;
}

vector<vector<int>> state_count(const vector<vector<uint8_t>> &data,
                                vector<int> &children, vector<int> &parents,
                                vector<int> &n_states, int &parallel) {
  // if parallel == 0 then use single thread, if parallel == 1 then use multi
  // thread(CPU), if parallel == 2 then use GPU (unimplemented)
  if (children.size() == 1) {
    int node_x = children.at(0);
    if (parents.empty()) {
      int x = n_states.at(node_x);
      vector<vector<int>> counts(x, vector<int>(1, 0));
      if (parallel == 0) {
        for (int i = 0; i < (int)data.size(); i++) {
          counts.at(data.at(i).at(node_x)).at(0) += 1;
        }
      } else if (parallel == 1) {
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
      } else if (parallel == 2) {
        // double* data_novec = data.data();
        // double* result = vec_res.data();
        // int n1 = vec_data.size();
        // int n2 = vec_data.at(0).size();

        // #pragma omp target data map(tofrom:data[0:n],result[0:n])

        // #pragma omp target data map(to:
        // data[0:data.size()][0:data.at(0).size()])
        // map(tofrom:counts[0:x][0:1]) #pragma omp target teams distribute
        // parallel
        // {
        //   vector<int> temp(x, 0);
        //   #pragma omp for thread_limit(256)
        //     for(int i = 0; i < data.size(); i++) {
        //       temp.at(data.at(i).at(node_x)) += 1;
        //     }
        //   for(int j = 0; j < x; j++){
        //     #pragma omp atomic
        //     counts.at(j).at(0) += temp.at(j);
        //   }
        // }
      }
      return counts;
    } else {
      // return the state counts of X, Y | Z shape: countmap[state of
      // child][state of parents]
      int x = n_states.at(node_x);
      int y = 1;
      for (int i = 0; i < (int)parents.size(); i++) {
        y = y * n_states.at(parents.at(i));
      }
      vector<vector<int>> counts(x, vector<int>(y, 0));
      // count the number of each state
      if (parallel == 0) {
        for (int i = 0; i < (int)data.size(); i++) {
          int yy = 0;
          for (int j = 0; j < (int)parents.size(); j++) {
            if (j == 0) {
              yy = data.at(i).at(parents.at(j));
            } else {
              yy = n_states.at(parents.at(j)) * yy +
                   data.at(i).at(parents.at(j));
            }
          }
          counts.at(data.at(i).at(node_x)).at(yy) += 1;
        }
      } else if (parallel == 1) {
        int yy;
#pragma omp parallel private(yy)
        {
          vector<vector<int>> temp(x, vector<int>(y, 0));
#pragma omp for
          for (int i = 0; i < (int)data.size(); i++) {
            yy = data.at(i).at(parents.at(0));
            for (int j = 1; j < (int)parents.size(); j++) {
              yy = n_states.at(parents.at(j)) * yy +
                   data.at(i).at(parents.at(j));
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
      } else if (parallel == 2) {
        return counts;
      }

      return counts;
    }
  } else {  // if children.size() == 2
    if (parents.empty()) {
      int xc = n_states.at(children.at(0));
      int yc = n_states.at(children.at(1));
      int len;
      len = xc * yc;
      vector<vector<int>> counts(len, vector<int>(1, 0));
      if (parallel == 0) {
        for (int i = 0; i < (int)data.size(); i++) {
          counts
              .at(data.at(i).at(children.at(0)) * yc +
                  data.at(i).at(children.at(1)))
              .at(0) += 1;
        }
      } else if (parallel == 1) {
#pragma omp parallel
        {
          vector<int> temp(len, 0);
#pragma omp for
          for (int i = 0; i < (int)data.size(); i++) {
            temp.at(data.at(i).at(children.at(0)) * yc +
                    data.at(i).at(children.at(1))) += 1;
          }
          for (int j = 0; j < len; j++) {
#pragma omp atomic
            counts.at(j).at(0) += temp.at(j);
          }
        }
      } else if (parallel == 2) {
        return counts;
      }

      return counts;
    } else {
      // return the state counts of X, Y | Z shape: countmap[state of
      // child][state of parents]
      int xc = n_states.at(children.at(0));
      int yc = n_states.at(children.at(1));
      int len;
      len = xc * yc;
      int y = 1;
      for (int i = 0; i < (int)parents.size(); i++) {
        y = y * n_states.at(parents.at(i));
      }
      vector<vector<int>> counts(len, vector<int>(y, 0));
      // count the number of each state
      if (parallel == 0) {
        for (int i = 0; i < (int)data.size(); i++) {
          int yy = 0;
          for (int j = 0; j < (int)parents.size(); j++) {
            if (j == 0) {
              yy = data.at(i).at(parents.at(j));
            } else {
              yy = n_states.at(parents.at(j)) * yy +
                   data.at(i).at(parents.at(j));
            }
          }
          counts
              .at(data.at(i).at(children.at(0)) * yc +
                  data.at(i).at(children.at(1)))
              .at(yy) += 1;
        }
      } else if (parallel == 1) {
        int yy;
#pragma omp parallel private(yy)
        {
          vector<vector<int>> temp(len, vector<int>(y, 0));
#pragma omp for
          for (int i = 0; i < (int)data.size(); i++) {
            yy = data.at(i).at(parents.at(0));
            for (int j = 1; j < (int)parents.size(); j++) {
              yy = n_states.at(parents.at(j)) * yy +
                   data.at(i).at(parents.at(j));
            }
            temp.at(data.at(i).at(children.at(0)) * yc +
                    data.at(i).at(children.at(1)))
                .at(yy) += 1;
          }
          for (int j = 0; j < len; j++) {
            for (int k = 0; k < y; k++) {
#pragma omp atomic
              counts.at(j).at(k) += temp.at(j).at(k);
            }
          }
        }
      } else if (parallel == 2) {
        return counts;
      }

      return counts;
    }
  }
}

vector<int> make_count_DP(const vector<vector<uint8_t>> &data, vector<int> &Gs,
                          vector<int> &n_states, int &parallel) {
  // Gs全ノードの頻度表を計算
  int x_len = 1;
  for (int i = 0; i < (int)Gs.size(); i++) {
    x_len = x_len * n_states.at(Gs.at(i));
  }
  vector<int> count_DP(x_len, 0);
  int yy;
#pragma omp parallel private(yy)
  {
    vector<int> temp(x_len, 0);
#pragma omp for
    for (int i = 0; i < (int)data.size(); i++) {
      yy = 0;
      for (int j = 0; j < (int)Gs.size(); j++) {
        if (j == 0) {
          yy = data.at(i).at(Gs.at(0));
        } else {
          yy = n_states.at(Gs.at(j)) * yy + data.at(i).at(Gs.at(j));
        }
      }
      temp.at(yy) += 1;
    }
    for (int j = 0; j < x_len; j++) {
#pragma omp atomic
      count_DP.at(j) += temp.at(j);
    }
  }
  return count_DP;
}

vector<vector<int>> state_count_DP(vector<int> &children, vector<int> &parents,
                                   vector<int> &n_states, vector<int> &count_DP,
                                   vector<int> &Gs, int &parallel) {
  // count the number of each state using dynamic programming
  int xx = 1;
  int yy = 1;
  for (int i = 0; i < (int)children.size(); i++) {
    xx = xx * n_states.at(children.at(i));
  }
  for (int i = 0; i < (int)parents.size(); i++) {
    yy = yy * n_states.at(parents.at(i));
  }
  vector<vector<int>> counts(xx, vector<int>(yy, 0));
  int size_Gs = Gs.size();
  vector<int> children_Gs;
  vector<int> parents_Gs;
  for (int i = 0; i < (int)children.size(); i++) {
    for (int j = 0; j < size_Gs; j++) {
      if (children.at(i) == Gs.at(j)) {
        children_Gs.push_back(j);
        break;
      }
    }
  }
  for (int i = 0; i < (int)parents.size(); i++) {
    for (int j = 0; j < size_Gs; j++) {
      if (parents.at(i) == Gs.at(j)) {
        parents_Gs.push_back(j);
        break;
      }
    }
  }
  vector<int> Gs_val(size_Gs, 0);
  if (parallel == 1) {
    int val, j, k, x_val, y_val;
#pragma omp parallel private(Gs_val, j, k, val, x_val, y_val)
    {
      Gs_val = vector<int>(size_Gs, 0);
      vector<vector<int>> temp(xx, vector<int>(yy, 0));
#pragma omp for
      for (int i = 0; i < (int)count_DP.size(); i++) {
        val = i;
        for (j = size_Gs - 1; j >= 0; j--) {
          Gs_val.at(j) = val % n_states.at(Gs.at(j));
          val = val / n_states.at(Gs.at(j));
        }  // calculate Gs_val from index of count_DP
        x_val = Gs_val.at(children_Gs.at(0));
        y_val = 0;

        for (j = 1; j < (int)children.size(); j++) {
          x_val = n_states.at(children.at(j)) * x_val +
                  Gs_val.at(children_Gs.at(j));
        }
        if (!parents.empty()) {
          y_val = Gs_val.at(parents_Gs.at(0));
        }
        for (j = 1; j < (int)parents.size(); j++) {
          y_val =
              n_states.at(parents.at(j)) * y_val + Gs_val.at(parents_Gs.at(j));
        }
        temp.at(x_val).at(y_val) += count_DP.at(i);
      }
      for (j = 0; j < xx; j++) {
        for (k = 0; k < yy; k++) {
#pragma omp atomic
          counts.at(j).at(k) += temp.at(j).at(k);
        }
      }
    }
  } else {
    for (int i = 0; i < (int)count_DP.size(); i++) {
      int val = i;
      for (int j = size_Gs - 1; j >= 0; j--) {
        Gs_val.at(j) = val % n_states.at(Gs.at(j));
        val = val / n_states.at(Gs.at(j));
      }  // calculate Gs_val from index of count_DP
      int x_val = Gs_val.at(children_Gs.at(0));
      int y_val = 0;

      for (int j = 1; j < (int)children.size(); j++) {
        x_val =
            n_states.at(children.at(j)) * x_val + Gs_val.at(children_Gs.at(j));
      }
      if (!parents.empty()) {
        y_val = Gs_val.at(parents_Gs.at(0));
      }
      for (int j = 1; j < (int)parents.size(); j++) {
        y_val =
            n_states.at(parents.at(j)) * y_val + Gs_val.at(parents_Gs.at(j));
      }
      counts.at(x_val).at(y_val) += count_DP.at(i);
    }
  }

  return counts;
}

float natori_independent_score(const vector<vector<uint8_t>> &data, int &node_x,
                               int &node_y, vector<int> &parents,
                               vector<int> &n_states, float &ESS, int &parallel,
                               bool &count_DP_flag, vector<int> &count_DP,
                               vector<int> &Gs) {
  // return log of the BDeu score of X, Y | Z
  double score = 0.0;
  double alpha = 0.5;  // hyperparameter
  vector<int> node_x_vec(1, node_x);
  vector<int> node_y_vec(1, node_y);
  if (parents.empty()) {
    // no parents
    vector<vector<int>> count;
    if (count_DP_flag) {
      count =
          state_count_DP(node_x_vec, parents, n_states, count_DP, Gs, parallel);
    } else {
      count = state_count(data, node_x_vec, parents, n_states, parallel);
    }
    int r = count.size();  // number of states of node_x
    int n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count.at(k).at(0);
    }
    for (int k = 0; k < r; k++) {  // for each state of node_x
      score += lgamma(alpha + count.at(k).at(0)) - lgamma(alpha);
    }
    score += lgamma(r * alpha) - lgamma(r * alpha + n_i);
    vector<vector<int>> count2;
    if (count_DP_flag) {
      count2 =
          state_count_DP(node_y_vec, parents, n_states, count_DP, Gs, parallel);
    } else {
      count2 = state_count(data, node_y_vec, parents, n_states, parallel);
    }
    r = count2.size();  // number of states of node_x
    n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count2.at(k).at(0);
    }
    for (int k = 0; k < r; k++) {  // for each state of node_x
      score += lgamma(alpha + count2.at(k).at(0)) - lgamma(alpha);
    }
    score += lgamma(r * alpha) - lgamma(r * alpha + n_i);
  } else {
    // have parents
    vector<vector<int>> count;
    if (count_DP_flag) {
      count =
          state_count_DP(node_x_vec, parents, n_states, count_DP, Gs, parallel);
    } else {
      count = state_count(data, node_x_vec, parents, n_states, parallel);
    }
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
        score += lgamma(alpha + count.at(k).at(j)) - lgamma(alpha);
        if (isnan(score)) {
          cout << ESS / (r * q) + count.at(k).at(j) << "score is nan" << endl;
        }
      }
      score += lgamma(r * alpha) - lgamma(r * alpha + n_ij.at(j));
    }
    vector<vector<int>> count2;
    if (count_DP_flag) {
      count2 =
          state_count_DP(node_y_vec, parents, n_states, count_DP, Gs, parallel);
    } else {
      count2 = state_count(data, node_y_vec, parents, n_states, parallel);
    }
    q = count2.at(0).size();  // number of states of parents
    r = count2.size();        // number of states of node_x
    vector<float> n_ij2(q, 0.0);
    for (int k = 0; k < r; k++) {
      for (int j = 0; j < q; j++) {
        n_ij2.at(j) += count2.at(k).at(j);
      }
    }
    for (int j = 0; j < q; j++) {    // for each state of parents
      for (int k = 0; k < r; k++) {  // for each state of node_x
        score += lgamma(alpha + count2.at(k).at(j)) - lgamma(alpha);
        if (isnan(score)) {
          cout << ESS / (r * q) + count2.at(k).at(j) << "score is nan" << endl;
        }
      }
      score += lgamma(r * alpha) - lgamma(r * alpha + n_ij2.at(j));
    }
  }
  // calculate the score
  return score;
}

float natori_dependent_score(const vector<vector<uint8_t>> &data, int &node_x,
                             int &node_y, vector<int> &parents,
                             vector<int> &n_states, float &ESS, int &parallel,
                             bool &count_DP_flag, vector<int> &count_DP,
                             vector<int> &Gs) {
  // return log of the BDeu score of X, Y | Z
  double score = 0.0;
  double alpha = 0.5;

  vector<int> children{node_x, node_y};
  if (parents.empty()) {
    // no parents
    vector<vector<int>> count;
    if (count_DP_flag) {
      count =
          state_count_DP(children, parents, n_states, count_DP, Gs, parallel);
    } else {
      count = state_count(data, children, parents, n_states, parallel);
    }
    int r = count.size();  // number of states of node_x]
    if (ESS > -1) {
      alpha = ESS / r;
    }
    int n_i = 0;
    for (int k = 0; k < r; k++) {
      n_i += count.at(k).at(0);
    }
    for (int k = 0; k < r; k++) {  // for each state of node_x
      score += lgamma(alpha + count.at(k).at(0)) - lgamma(alpha);
    }
    score += lgamma(r * alpha) - lgamma(r * alpha + n_i);
  } else {
    // have parents
    vector<vector<int>> count;
    if (count_DP_flag) {
      count =
          state_count_DP(children, parents, n_states, count_DP, Gs, parallel);
    } else {
      count = state_count(data, children, parents, n_states, parallel);
    }
    int q = count.at(0).size();  // number of states of parents
    int r = count.size();        // number of states of node_x
    if (ESS > -1) {
      alpha = ESS / (r * q);
    }
    vector<float> n_ij(q, 0.0);
    for (int k = 0; k < r; k++) {
      for (int j = 0; j < q; j++) {
        n_ij.at(j) += count.at(k).at(j);
      }
    }
    for (int j = 0; j < q; j++) {    // for each state of parents
      for (int k = 0; k < r; k++) {  // for each state of node_x
        score += lgamma(alpha + count.at(k).at(j)) - lgamma(alpha);
        if (isnan(score)) {
          cout << ESS / (r * q) + count.at(k).at(j) << "score is nan" << endl;
        }
      }
      score += lgamma(r * alpha) - lgamma(r * alpha + n_ij.at(j));
    }
  }
  return score;
}

bool ci_test(const vector<vector<uint8_t>> &data, int &node_x, int &node_y,
             vector<int> &Z, vector<int> &n_states, float &ESS, int &parallel,
             bool &count_DP_flag, vector<int> &count_DP, vector<int> &Gs,
             int &n_citest, int &n_citest_DP, vector<vector<vector<bool>>> &sep,
             vector<vector<bool>> &sep2) {
  // CI test for X _|_ Y | Z
  if (sep2.at(node_x).at(node_y)) {
    return true;
  }
  float independent_score = 0.0;
  float dependent_score = 0.0;
  independent_score +=
      natori_independent_score(data, node_x, node_y, Z, n_states, ESS, parallel,
                               count_DP_flag, count_DP, Gs);
  dependent_score +=
      natori_dependent_score(data, node_x, node_y, Z, n_states, ESS, parallel,
                             count_DP_flag, count_DP, Gs);

  if (count_DP_flag) {
    n_citest_DP += 1;
  } else {
    n_citest += 1;
  }
  if (independent_score > dependent_score) {
    sep2.at(node_x).at(node_y) = true;
    sep2.at(node_y).at(node_x) = true;
    for (int j = 0; j < (int)Z.size(); j++) {
      sep.at(node_x).at(node_y).at(Z.at(j)) = true;
      sep.at(node_y).at(node_x).at(Z.at(j)) = true;
    }
    // cout<< "CI independent:" <<node_x<<" _|_"<<node_y<<"
    // |"<<independent_score<<">"<<dependent_score<< endl;
    return true;
  } else {
    // cout<< "CI dependent:" <<node_x<<" _|_"<<node_y<<"
    // |"<<independent_score<<"<"<<dependent_score<< endl;
    return false;
  }
}

void orientation_A2(PDAG_nishikori &Gall, vector<int> &Gs,
                    vector<vector<bool>> &deletededges) {
  /*
      orient edges in a PDAG_nishikori to a maximally oriented graph.
      orient rules are based on rule 1~3 from Meek,C.:Causal Inference and
     Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in
     Artificial Inteligence (UAl-95),p.403-410 (195) in this stage (stage A2 in
     Yehezkel and Lerner(2009)), only rule 1 is applied because only X -> Y - Z
     shape is created in stage A1 (X-Z removed).
  */
  // Rule 1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
  for (int i = 0; i < (int)Gall.g.size(); i++) {
    for (int j = 0; j < (int)Gall.g.size(); j++) {
      if (deletededges.at(i).at(j) || deletededges.at(j).at(i)) {
        int X = i;
        int Z = j;
        for (auto &Y : Gall.undirected_neighbors(Z)) {
          if (Gall.has_directed_edge(X, Y)) {
            Gall.remove_edge(Z, Y);
            // cout<< "A2_oriented:" <<Y<<"->"<<Z<< endl;
          }
        }
      }
    }
  }
  return;
}

void orientation_B2(PDAG_nishikori &Gall, vector<int> &Gs,
                    vector<vector<bool>> &deletededges,
                    const vector<vector<uint8_t>> &data, vector<int> &n_states,
                    float &ESS, int &parallel, bool &count_DP_flag,
                    vector<int> &count_DP, vector<vector<vector<bool>>> &sep,
                    vector<vector<bool>> &sep2) {
  /*
      orient edges in a PDAG_nishikori to a maximally oriented graph.
      orient rules are based on rule 1~3 from Meek,C.:Causal Inference and
     Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in
     Artificial Inteligence (UAl-95),p.403-410 (195)
  */
  // for each X-Z-Y (X and Y is not adjecent), find V-structure and orient as X
  // -> Z <- Y
  bool store = count_DP_flag;
  vector<bool> Gs_flag(Gall.g.size(), true);
  for (auto &node : Gs) {
    Gs_flag.at(node) = true;
  }
  vector<vector<bool>> vstructuredetected(Gall.g.size(),
                                          vector<bool>(Gall.g.size(), false));
  for (auto &X : Gs) {
    for (auto &Z : Gall.undirected_neighbors(X)) {
      for (auto &Y : Gall.undirected_neighbors(Z)) {
        // if (X != Y && !Gall.has_edge(X, Y) && !Gall.has_edge(Y, X) &&
        // Gall.has_edge(X, Z) && Gall.has_edge(Z, X) && Gall.has_edge(Y, Z) &&
        // Gall.has_edge(Z, Y)) {
        if (X != Y && !Gall.has_edge(X, Y) && !Gall.has_edge(Y, X) &&
            (!vstructuredetected.at(X).at(Z) ||
             !vstructuredetected.at(Y).at(Z))) {
          if (Gs_flag.at(Y) == false || Gs_flag.at(Z) == false) {
            count_DP_flag = false;
          }
          // cout<< "V-structure think:" <<X<<"->"<<Z<<"<-"<<Y<< endl;
          if (!sep.at(X).at(Y).at(Z) && sep2.at(X).at(Y)) {
            vstructuredetected.at(X).at(Z) = true;
            vstructuredetected.at(Y).at(Z) = true;
            // cout<< "V-structure found:" <<X<<"->"<<Z<<"<-"<<Y<< endl;
          }
          count_DP_flag = store;
        }
      }
    }
  }
  for (int i = 0; i < (int)Gall.g.size(); i++) {
    for (int j = 0; j < (int)Gall.g.size(); j++) {
      if (vstructuredetected.at(i).at(j) && !vstructuredetected.at(j).at(i)) {
        Gall.remove_edge(j, i);
      }
    }
  }
  count_DP_flag = store;
  bool flag = true;
  while (flag) {
    flag = false;
    // Rule 1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
    for (auto &X : Gs) {
      for (auto &Y : Gall.successors(X)) {
        if (!Gall.has_edge(Y, X) && Gall.has_edge(X, Y)) {
          for (auto &Z : Gall.undirected_neighbors(Y)) {
            if (!Gall.has_edge(X, Z) && !Gall.has_edge(Z, X) && Z != X) {
              Gall.remove_edge(Z, Y);
              // cout<< "R1:" <<Y<<"->"<<Z<<"|"<<X<< endl;
              flag = true;
            }
          }
        }
      }
    }
    // // Rule 1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
    // for (auto &X : Gs) {
    //   for (auto &Y : Gs) {
    //     for (auto &Z : Gs) {
    //       if (!Gall.has_edge(X, Z) && !Gall.has_edge(Z, X) && Z != X && Z !=
    //       Y && X != Y && Gall.has_edge(X, Y) && !Gall.has_edge(Y, X) &&
    //       Gall.has_edge(Y, Z) && Gall.has_edge(Z, Y)) {
    //         Gall.remove_edge(Z, Y);
    //         // cout<< "R1:" <<Y<<"->"<<Z<<"|"<<X<< endl;
    //         flag = true;
    //       }
    //     }
    //   }
    // }
    // Rule 2: X - Y and if there is a directed path from X to Y, then X -> Y
    for (auto &X : Gs) {
      for (auto &Y : Gall.undirected_neighbors(X)) {
        if (Gall.has_directed_path(X, Y)) {
          Gall.remove_edge(Y, X);
          // cout<< "R2:" <<X<<"->"<<Y<< endl;
          flag = true;
        }
      }
    }
    // Rule 3: for each X->W<-Z X-Y-Z Y-W, orient Y->W
    for (auto &X : Gs) {
      for (auto &Y : Gall.undirected_neighbors(X)) {
        for (auto &Z : Gall.undirected_neighbors(Y)) {
          if (Z != X && !Gall.has_edge(Z, X) &&
              !Gall.has_edge(X, Z)) {  // X-Y-Z
            for (auto &W : Gall.undirected_neighbors(Y)) {
              if (W != X && W != Z && Gall.has_directed_edge(X, W) &&
                  Gall.has_directed_edge(Z, W)) {
                Gall.remove_edge(W, Y);
                // cout<< "R3:" <<Y<<"->"<<W<<"|"<<"X="<<X<<"Z="<<Z<< endl;
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

void recursive_search(const vector<vector<uint8_t>> &data, PDAG_nishikori &Gall,
                      vector<int> Gs, vector<int> &Gex, int N,
                      vector<int> &n_states, float &ESS, int &parallel,
                      int &threshold_DP, bool &search_neighbor,
                      bool &do_orientation_A2,
                      vector<vector<vector<bool>>> &sep,
                      vector<vector<bool>> &sep2, int &n_citest,
                      int &n_citest_DP, int &n_DP) {
  int n_node = data.at(0).size();
  vector<int> count_DP;
  bool count_DP_flag = false;

  // stage 0: exit condition
  bool exitcondition = true;
  for (int i = 0; i < (int)Gs.size(); i++) {
    if ((int)Gall.predecessors(Gs.at(i)).size() > N) {
      exitcondition = false;
      break;
    }
  }
  if (exitcondition) {
    return;
  }
  // stage A1: Do CI tests for nodes between Gex and Gs and remove edge
  vector<vector<bool>> deletededges(n_node, vector<bool>(n_node, false));
  if (!(Gex.empty())) {
    for (auto &node_y : Gs) {
      for (auto &node_x : Gex) {
        if (Gall.has_edge(node_x, node_y) &&
            !(deletededges.at(node_x).at(node_y) ||
              deletededges.at(node_y).at(node_x))) {
          vector<int> Z = Gall.predecessors(node_y);
          for (int i = 0; i < (int)Z.size(); i++) {
            if (Z.at(i) == node_x) {
              Z.erase(Z.begin() + i);
            }
          }  // erase node_x from Z
          sort(Z.begin(), Z.end());
          if ((int)Z.size() >= N) {
            do {
              vector<int> selected_z;
              for (int j = 0; j < N; j++) {
                selected_z.push_back(Z.at(j));
              }
              if (!deletededges.at(node_x).at(node_y)) {
                if (ci_test(data, node_x, node_y, selected_z, n_states, ESS,
                            parallel, count_DP_flag, count_DP, Gs, n_citest,
                            n_citest_DP, sep, sep2)) {
                  Gall.remove_edge(node_x, node_y);
                  deletededges.at(node_x).at(node_y) = true;
                  deletededges.at(node_y).at(node_x) = true;
                  // cout<< "A1_removed:" <<node_x<<"-"<<node_y<< endl;
                  // transive_cut();
                }
              }
            } while (next_combination(Z.begin(), Z.end(), N));
          }
        }
      }
    }
  }
  // stage A2: orient edges in Gs using "smart" orientation rules R1
  if (do_orientation_A2) {
    orientation_A2(Gall, Gs, deletededges);
  }

  // stage B1: Do CI tests for nodes between Gs and Gs and remove edge

  // make DP map if the number of nodes in GS is smaller than threshold_DP

  if (threshold_DP > (int)Gs.size() && threshold_DP != 0 && Gs.size() > 2) {
    count_DP = make_count_DP(data, Gs, n_states, parallel);
    n_DP += 1;
    count_DP_flag = true;
  }
  vector<bool> Gs_flag(n_node, false);
  for (auto &node : Gs) {
    Gs_flag.at(node) = true;
  }
  for (auto &node_y : Gs) {
    for (auto &node_x : Gs) {
      if (!deletededges.at(node_x).at(node_y) && node_x != node_y &&
          (Gall.has_edge(node_x, node_y) || Gall.has_edge(node_y, node_x))) {
        if (N == 0) {
          vector<int> S;
          if (ci_test(data, node_x, node_y, S, n_states, ESS, parallel,
                      count_DP_flag, count_DP, Gs, n_citest, n_citest_DP, sep,
                      sep2)) {
            Gall.remove_edge(node_x, node_y);
            Gall.remove_edge(node_y, node_x);
            deletededges.at(node_x).at(node_y) = true;
            deletededges.at(node_y).at(node_x) = true;
            // cout<< "B1_removed:" <<node_x<<"-"<<node_y<< endl;
            // transive_cut();
          }
        } else {
          vector<int> S;
          if (search_neighbor) {
            S = Gall.neighbors(node_y);
          } else {
            S = Gall.predecessors(node_y);
          }
          auto newEnd = remove(S.begin(), S.end(), node_x);
          S.erase(newEnd, S.end());  // remove node_x from S  正しい？
          if ((int)S.size() >= N) {
            do {
              vector<int> selected_z;
              bool store_flag = count_DP_flag;
              for (int j = 0; j < N; j++) {
                selected_z.push_back(S.at(j));
                if (!Gs_flag.at(selected_z.at(j))) {
                  count_DP_flag =
                      false;  // check if DP can be used for this CI test
                }
              }
              if (ci_test(data, node_x, node_y, selected_z, n_states, ESS,
                          parallel, count_DP_flag, count_DP, Gs, n_citest,
                          n_citest_DP, sep, sep2)) {
                Gall.remove_edge(node_x, node_y);
                Gall.remove_edge(node_y, node_x);
                deletededges.at(node_x).at(node_y) = true;
                deletededges.at(node_y).at(node_x) = true;
                // cout<< "B1_removed:" <<node_x<<"-"<<node_y<< endl;
                // transive_cut();
              }
              count_DP_flag = store_flag;
            } while (next_combination(S.begin(), S.end(), N));
          }
        }
      }
    }
  }

  // stage B2: orient edges in Gs using orientation rules R1~R3
  orientation_B2(Gall, Gs, deletededges, data, n_states, ESS, parallel,
                 count_DP_flag, count_DP, sep, sep2);
  // stage B3: Group the nodes having the lowest topological order into a
  // descendant substructure Gd
  vector<int> Gd;
  vector<vector<int>> g_ex_connection;
  order_grouping(Gall, Gs, Gd, g_ex_connection);

  // stage C: Ancestor sub-structure decomposition
  vector<int> Gexd;
  for (int i = 0; i < (int)g_ex_connection.size(); i++) {
    for (int j = 0; j < (int)g_ex_connection.at(i).size(); j++) {
      Gexd.push_back(g_ex_connection.at(i).at(j));
    }
    recursive_search(data, Gall, g_ex_connection.at(i), Gex, N + 1, n_states,
                     ESS, parallel, threshold_DP, search_neighbor,
                     do_orientation_A2, sep, sep2, n_citest, n_citest_DP, n_DP);
  }
  sort(Gexd.begin(), Gexd.end());
  for (int i = 0; i < (int)Gexd.size(); i++) {
    bool check = true;
    for (int j = 0; j < (int)Gex.size(); j++) {
      if (Gexd.at(i) == Gex.at(j)) {
        check = false;
        break;
      }
    }
    if (check) {
      Gex.push_back(Gexd.at(i));
    }
  }

  // stage D:  Descendantsub-structure decomposition
  recursive_search(data, Gall, Gd, Gex, N + 1, n_states, ESS, parallel,
                   threshold_DP, search_neighbor, do_orientation_A2, sep, sep2,
                   n_citest, n_citest_DP, n_DP);
  return;
}

py::array_t<bool> RAI_nishikori(py::array_t<uint8_t> data,
                                py::array_t<int> n_states,  // uint8_t
                                float ESS, int parallel, int threshold_DP,
                                bool search_neighbor, bool do_orientation_A2) {
  // translate imput data to c++ vector(this is not optimal but I don't know how
  // to use pybind11::array_t)
  py::buffer_info buf_data = data.request(), buf_states = n_states.request();
  const uint8_t *__restrict__ prt_data = static_cast<uint8_t *>(buf_data.ptr);
  const int *__restrict__ prt_states = static_cast<int *>(buf_states.ptr);
  int n_data = buf_data.shape[0],
      n_node = buf_data.shape[1];  // number of nodes
  vector<vector<uint8_t>> data_vec(n_data, vector<uint8_t>(n_node));
  vector<int> n_states_vec(n_node, 0);
  for (int i = 0; i < n_data; i++) {
    for (int j = 0; j < n_node; j++) {
      data_vec.at(i).at(j) = prt_data[i * n_node + j];
    }
  }
  for (int i = 0; i < n_node; i++) {
    n_states_vec.at(i) = prt_states[i];
  }
  auto endg = py::array_t<bool>({n_node, n_node});
  py::buffer_info buf_endg = endg.request();
  bool *__restrict__ prt_endg = static_cast<bool *>(buf_endg.ptr);

  int n_citest = 0;
  int n_citest_DP = 0;
  int n_DP = 0;
  int sep_size = n_node;
  vector<vector<vector<bool>>> sep(
      sep_size, vector<vector<bool>>(sep_size, vector<bool>(sep_size, false)));
  vector<vector<bool>> sep2(sep_size, vector<bool>(sep_size, false));

  // initialize Gall, Gs, Gex
  PDAG_nishikori Gall;
  vector<vector<bool>> gall(n_node, vector<bool>(n_node, true));
  for (int i = 0; i < n_node; i++) {
    gall.at(i).at(i) = false;
  }
  Gall.g = gall;  // Gall is complete graph
  vector<int> Gs(n_node, -1);
  for (int i = 0; i < n_node; i++) {
    Gs.at(i) = i;
  }                 // Gs contains all nodes 0 ~ n_node - 1
  vector<int> Gex;  // Gex is empty
  PDAG_nishikori Gend;

  recursive_search(data_vec, Gall, Gs, Gex, 0, n_states_vec, ESS, parallel,
                   threshold_DP, search_neighbor, do_orientation_A2, sep, sep2,
                   n_citest, n_citest_DP, n_DP);

  // translate Gall to py::array_t (this is not optimal but I don't know how to
  // use pybind11::array_t)
  for (int i = 0; i < n_node; i++) {
    for (int j = 0; j < n_node; j++) {
      prt_endg[i * n_node + j] = Gall.g.at(i).at(j);
    }
  }
  cout << "citest_count: " << n_citest << ", citest_DP_count: " << n_citest_DP
       << ", DPmap_count: " << n_DP << endl;
  return endg;
}