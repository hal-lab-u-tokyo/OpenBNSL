#include <cassert>
#include <cstddef>
#include <set>
#include <stack>
#include <string>
#include <vector>
using namespace std;

#include "structure_learning/cuPC.h"
// the following includes are for permutation and combination algorithms
#include <algorithm>
#include <functional>
#include <iostream>

// for gamma function
#include <cmath>

// input: data: np.ndarray,  shape: (n: number of variables, d: number of
// samples) output: leard PDAG
// Gall.at(i).at(j)==1 means there is an edge i -> j

namespace cuda3 {
#define CUDA_CHECK(call)                               \
  do {                                                 \
    cudaError_t e = call;                              \
    if (e != cudaSuccess) {                            \
      throw std::runtime_error(cudaGetErrorString(e)); \
    }                                                  \
  } while (0)

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

const int max_level = 20;
const int max_dim = 21;

// Behrooz Zarebavani, Foad Jafarinejad, Matin Hashemi, Saber Salehkaleybar,
// cuPC: CUDA-based Parallel PC Algorithm for Causal Structure Learning on GPU,
// IEEE Transactions on Parallel and Distributed Systems (TPDS), Vol. 31, No. 3,
// March 2020.

__device__ double PYTHAG(double a, double b) {
  double at = abs(a), bt = abs(b), ct, result;

  if (at > bt) {
    ct = bt / at;
    result = at * sqrt(1.0 + ct * ct);
  } else if (bt > 0.0) {
    ct = at / bt;
    result = bt * sqrt(1.0 + ct * ct);
  } else {
    result = 0.0;
  }
  return (result);
}

__device__ double SIGN(double a, double b) {
  return b >= 0.0 ? abs(a) : -abs(a);
}

__device__ void pseudoinverse(int m, double M2[][max_level],
                              double M2Inv[][max_level]) {
  double v[max_level][max_level];
  double rv1[max_level];
  double w[max_level];
  double res1[max_level][max_level];
  int flag, its, i, j, jj, k, l, nm;
  double c, f, h, s, x, y, z;
  double anorm = 0.0, g = 0.0, scale = 0.0;
  /* Householder reduction to bidiagonal form */
  for (i = 0; i < m; i++) {
    /* left-hand reduction */
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;
    if (i < m) {
      for (k = i; k < m; k++) scale += abs(M2[k][i]);
      if (scale) {
        for (k = i; k < m; k++) {
          M2[k][i] = (M2[k][i] / scale);
          s += (M2[k][i] * M2[k][i]);
        }
        f = M2[i][i];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        M2[i][i] = f - g;
        if (i != m - 1) {
          for (j = l; j < m; j++) {
            for (s = 0.0, k = i; k < m; k++) s += (M2[k][i] * M2[k][j]);
            f = s / h;
            for (k = i; k < m; k++) M2[k][j] += (f * M2[k][i]);
          }
        }
        for (k = i; k < m; k++) M2[k][i] = (M2[k][i] * scale);
      }
    }
    w[i] = scale * g;

    /* right-hand reduction */
    g = s = scale = 0.0;
    if (i < m && i != m - 1) {
      for (k = l; k < m; k++) scale += abs(M2[i][k]);
      if (scale) {
        for (k = l; k < m; k++) {
          M2[i][k] = (M2[i][k] / scale);
          s += (M2[i][k] * M2[i][k]);
        }
        f = M2[i][l];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        M2[i][l] = f - g;
        for (k = l; k < m; k++) rv1[k] = M2[i][k] / h;
        if (i != m - 1) {
          for (j = l; j < m; j++) {
            for (s = 0.0, k = l; k < m; k++) s += (M2[j][k] * M2[i][k]);
            for (k = l; k < m; k++) M2[j][k] += (s * rv1[k]);
          }
        }
        for (k = l; k < m; k++) M2[i][k] = M2[i][k] * scale;
      }
    }
    anorm = max(anorm, (abs(w[i]) + abs(rv1[i])));
  }

  /* accumulate the right-hand transformation */
  for (i = m - 1; i >= 0; i--) {
    if (i < m - 1) {
      if (g) {
        for (j = l; j < m; j++) v[j][i] = (M2[i][j] / M2[i][l]) / g;
        /* double division to avoid underflow */
        for (j = l; j < m; j++) {
          for (s = 0.0, k = l; k < m; k++) s += (M2[i][k] * v[k][j]);
          for (k = l; k < m; k++) v[k][j] += (s * v[k][i]);
        }
      }
      for (j = l; j < m; j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }

  /* accumulate the left-hand transformation */
  for (i = m - 1; i >= 0; i--) {
    l = i + 1;
    g = w[i];
    if (i < m - 1)
      for (j = l; j < m; j++) M2[i][j] = 0.0;
    if (g) {
      g = 1.0 / g;
      if (i != m - 1) {
        for (j = l; j < m; j++) {
          for (s = 0.0, k = l; k < m; k++) s += (M2[k][i] * M2[k][j]);
          f = (s / M2[i][i]) * g;
          for (k = i; k < m; k++) M2[k][j] += (f * M2[k][i]);
        }
      }
      for (j = i; j < m; j++) M2[j][i] = (M2[j][i] * g);
    } else {
      for (j = i; j < m; j++) M2[j][i] = 0.0;
    }
    ++M2[i][i];
  }

  /* diagonalize the bidiagonal form */
  for (k = m - 1; k >= 0; k--) {     /* loop over singular values */
    for (its = 0; its < 30; its++) { /* loop over allowed iterations */
      flag = 1;
      for (l = k; l >= 0; l--) { /* test for splitting */
        nm = l - 1;
        if (abs(rv1[l]) + anorm == anorm) {
          flag = 0;
          break;
        }
        if (abs(w[nm]) + anorm == anorm) break;
      }
      if (flag) {
        c = 0.0;
        s = 1.0;
        for (i = l; i <= k; i++) {
          f = s * rv1[i];
          if (abs(f) + anorm != anorm) {
            g = w[i];
            h = PYTHAG(f, g);
            w[i] = h;
            h = 1.0 / h;
            c = g * h;
            s = (-f * h);
            for (j = 0; j < m; j++) {
              y = M2[j][nm];
              z = M2[j][i];
              M2[j][nm] = (y * c + z * s);
              M2[j][i] = (z * c - y * s);
            }
          }
        }
      }
      z = w[k];
      if (l == k) {    /* convergence */
        if (z < 0.0) { /* make singular value nonnegative */
          w[k] = (-z);
          for (j = 0; j < m; j++) v[j][k] = (-v[j][k]);
        }
        break;
      }
      if (its >= 30) {
        printf("Not converged\n");
      }

      /* shift from bottom 2 x 2 minor */
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = PYTHAG(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

      /* next QR transformation */
      c = s = 1.0;
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = PYTHAG(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;
        for (jj = 0; jj < m; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = (x * c + z * s);
          v[jj][i] = (z * c - x * s);
        }
        z = PYTHAG(f, h);
        w[j] = z;
        if (z) {
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }
        f = (c * g) + (s * y);
        x = (c * y) - (s * g);
        for (jj = 0; jj < m; jj++) {
          y = M2[jj][j];
          z = M2[jj][i];
          M2[jj][j] = (y * c + z * s);
          M2[jj][i] = (z * c - y * s);
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }

  // start compute inverse matrix

  for (int rowNumber = 0; rowNumber < m; rowNumber++) {
    for (int colNumber = 0; colNumber < m; colNumber++) {
      res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
    }
  }

  for (int rowNumber = 0; rowNumber < m; rowNumber++) {
    for (int colNumber = 0; colNumber < m; colNumber++) {
      M2Inv[rowNumber][colNumber] = 0;
      for (int thirdIndex = 0; thirdIndex < m; thirdIndex++) {
        M2Inv[rowNumber][colNumber] =
            M2Inv[rowNumber][colNumber] +
            res1[rowNumber][thirdIndex] * M2[colNumber][thirdIndex];
      }
    }
  }
}

__device__ int binom(int n, int k) {
  if (k > n - k) k = n - k;
  int res = 1;
  for (int i = 0; i < k; i++) {
    res *= n - i;
    res /= i + 1;
  }
  return res;
}

__device__ void comb(int n, int l, int t, int p, int *idset) {
  int sum = 0;
  int max_t = binom(n, l) - 1;
  if (t > max_t) t = max_t;
  for (int i = 0; i < l; i++) {
    int x = (i == 0 ? 0 : idset[i - 1] + 1);
    while (true) {
      int d = binom(n - x - 1, l - i - 1);
      if (sum + d > t) break;
      x++;
      sum += d;
    }
    idset[i] = x;
  }
  if (p >= 0) {
    for (int i = 0; i < l; i++) {
      if (idset[i] >= p) {
        idset[i]++;
      }
    }
  }
}

__global__ void calc_correlation_matrix(int n_node, int n_data, uint8_t *data,
                                        double *C) {
  int i = blockIdx.x;
  int j = blockIdx.y;
  double mean_i = 0, mean_j = 0, mean_i2 = 0, mean_j2 = 0, mean_ij = 0;
  for (int k = 0; k < n_data; k++) {
    double data_i = data[i * n_data + k];
    double data_j = data[j * n_data + k];
    mean_i += data_i;
    mean_j += data_j;
    mean_i2 += data_i * data_i;
    mean_j2 += data_j * data_j;
    mean_ij += data_i * data_j;
  }
  mean_i /= n_data;
  mean_j /= n_data;
  mean_i2 /= n_data;
  mean_j2 /= n_data;
  mean_ij /= n_data;
  double sigma_i = sqrt(mean_i2 - mean_i * mean_i);
  double sigma_j = sqrt(mean_j2 - mean_j * mean_j);
  if (sigma_i == 0 || sigma_j == 0) {
    C[i * n_node + j] = 0;
  } else {
    C[i * n_node + j] = (mean_ij - mean_i * mean_j) / (sigma_i * sigma_j);
  }
}

__global__ void PC_level_0(int n_node, int n_data, int *G, double *C) {
  int i = blockIdx.x;
  int j = blockIdx.y;
  if (i >= j) {
    return;
  }
  double rho = C[i * n_node + j];
  double z = abs(0.5 * log((1 + rho) / (1 - rho)));
  // double tau = 1.96 / sqrt(static_cast<double>(n_data - 3));
  double tau = 4 / sqrt(static_cast<double>(n_data - 3));
  if (z <= tau) {
    G[i * n_node + j] = 0;
    G[j * n_node + i] = 0;
  }
}

__global__ void PC_level_n(int level, int n_node, int n_data, int *G, double *C,
                           int *sepsets) {
  extern __shared__ int smem[];
  int i = blockIdx.x;
  int *G_compacted = smem;
  if (threadIdx.x == 0) {
    int cnt = 0;
    for (int j = 0; j < n_node; j++) {
      if (G[i * n_node + j] == 1) {
        G_compacted[++cnt] = j;
      }
    }
    G_compacted[0] = cnt;
  }
  __syncthreads();
  int n_adj = G_compacted[0];
  if (n_adj - 1 < level) {
    return;
  }
  int sepset_cnt = binom(n_adj - 1, level);
  int step = gridDim.y * blockDim.x;
  int sepset_cnt_loop = (sepset_cnt + step - 1) / step * step;
  int sepset[max_level];
  double H[2][2], M1[2][max_level], M2[max_level][max_level],
      M2Inv[max_level][max_level], M1mulM2Inv[2][max_level];
  for (int sepset_idx = blockIdx.y * blockDim.x + threadIdx.x;
       sepset_idx < sepset_cnt_loop; sepset_idx += step) {
    comb(n_adj, level, sepset_idx, -1, sepset);
    for (int k = 0; k < level; k++) {
      sepset[k] = G_compacted[sepset[k] + 1];
    }
    for (int k = 0; k < level; k++) {
      for (int l = 0; l < level; l++) {
        M2[k][l] = C[sepset[k] * n_node + sepset[l]];
      }
    }
    pseudoinverse(level, M2, M2Inv);
    for (int idx_j = 0; idx_j < n_adj; idx_j++) {
      int j = G_compacted[idx_j + 1];
      if (!G[i * n_node + j]) continue;
      bool in_sepset = false;
      for (int k = 0; k < level; k++) {
        if (j == sepset[k]) {
          in_sepset = true;
          break;
        }
      }
      if (in_sepset) {
        continue;
      }
      for (int k = 0; k < level; k++) {
        M1[0][k] = C[i * n_node + sepset[k]];
        M1[1][k] = C[j * n_node + sepset[k]];
      }
      for (int k = 0; k < 2; k++) {
        for (int l = 0; l < level; l++) {
          M1mulM2Inv[k][l] = 0;
          for (int m = 0; m < level; m++) {
            M1mulM2Inv[k][l] += M1[k][m] * M2Inv[m][l];
          }
        }
      }
      H[0][0] = H[1][1] = 1;
      H[0][1] = H[1][0] = C[i * n_node + j];
      for (int k = 0; k < 2; k++) {
        for (int l = 0; l < 2; l++) {
          for (int m = 0; m < level; m++) {
            H[k][l] -= M1mulM2Inv[k][m] * M1[l][m];
          }
        }
      }
      double rho = H[0][1] / sqrt(H[0][0] * H[1][1]);
      double z = abs(0.5 * log((1 + rho) / (1 - rho)));
      // double tau = 1.96 / sqrt(static_cast<double>(n_data - level - 3));
      double tau = 4 / sqrt(static_cast<double>(n_data - level - 3));
      if (z <= tau) {
        int ij_min = (i < j ? i : j);
        int ij_max = (i < j ? j : i);
        if (atomicCAS(G + ij_min * n_node + ij_max, 1, 0) == 1) {
          G[ij_max * n_node + ij_min] = 0;
          for (int k = 0; k < level; k++) {
            sepsets[(ij_min * n_node + ij_max) * max_level + k] = sepset[k];
          }
        }
      }
    }
  }
}

void orientation(PDAG &G, const vector<int> &sepsets) {
  /*
      orient edges in a PDAG to a maximally oriented graph.
      orient rules are based on rule 1~3 from Meek,C.:Causal Inference and
     Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in
     Artificial Inteligence (UAl-95),p.403-410 (195)
  */
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

PDAG PCsearch(int n_node, int n_data, const vector<uint8_t> &data,
              const vector<int> &n_states) {
  vector<int> G(n_node * n_node);
  for (int i = 0; i < n_node; i++) {
    for (int j = 0; j < n_node; j++) {
      if (i != j) G[i * n_node + j] = 1;
    }
  }
  vector<double> C(n_node * n_node);
  vector<int> sepsets(n_node * n_node * max_level, -1);
  int *G_d, *sepsets_d;
  uint8_t *data_d;
  double *C_d;
  int size_G = sizeof(int) * n_node * n_node;
  int size_data = sizeof(uint8_t) * n_data * n_node;
  int size_C = sizeof(double) * n_node * n_node;
  int size_sepsets = sizeof(int) * n_node * n_node * max_level;
  CUDA_CHECK(cudaMalloc(&G_d, size_G));
  CUDA_CHECK(cudaMalloc(&data_d, size_data));
  CUDA_CHECK(cudaMalloc(&C_d, size_C));
  CUDA_CHECK(cudaMalloc(&sepsets_d, size_sepsets));
  CUDA_CHECK(cudaMemcpy(G_d, G.data(), size_G, cudaMemcpyHostToDevice));
  CUDA_CHECK(
      cudaMemcpy(data_d, data.data(), size_data, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(C_d, C.data(), size_C, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(sepsets_d, sepsets.data(), size_sepsets,
                        cudaMemcpyHostToDevice));
  calc_correlation_matrix<<<dim3(n_node, n_node), 1>>>(n_node, n_data, data_d,
                                                       C_d);
  CUDA_CHECK(cudaMemcpy(C.data(), C_d, size_C, cudaMemcpyDeviceToHost));
  // cout << "correlation matrix:" << endl;
  // for (int i = 0; i < n_node; i++) {
  //   for (int j = 0; j < n_node; j++) {
  //     cout << C[i * n_node + j] << ' ';
  //   }
  //   cout << endl;
  // }
  // stage 1: Do CI tests between nodes and remove edge (undirected graph)
  int level = 0;
  int max_n_adj = n_node - 1;
  uint64_t max_dim_s = 1;
  while (level <= n_node - 2) {
    cout << "level: " << level << ", max_n_adj: " << max_n_adj << endl;
    if (level == 0) {
      dim3 threadsPerBlock(1);
      dim3 numBlocks(n_node, n_node);
      PC_level_0<<<numBlocks, threadsPerBlock>>>(n_node, n_data, G_d, C_d);
    } else {
      dim3 threadsPerBlock(64);
      dim3 numBlocks(n_node, 2);
      PC_level_n<<<numBlocks, threadsPerBlock, sizeof(int) * max_n_adj>>>(
          level, n_node, n_data, G_d, C_d, sepsets_d);
    }
    CUDA_CHECK(cudaMemcpy(G.data(), G_d, size_G, cudaMemcpyDeviceToHost));
    max_n_adj = 0;
    int next_node_cnt = 0;
    int next_edge_cnt = 0;
    for (int i = 0; i < n_node; i++) {
      int n_adj = 0;
      for (int j = 0; j < n_node; j++) {
        if (G[i * n_node + j]) n_adj++;
      }
      if (n_adj - 1 > level) {
        next_node_cnt++;
      }
      max_n_adj = max(max_n_adj, n_adj);
      next_edge_cnt += n_adj;
    }
    cout << "next_node_cnt, next_edge_cnt: " << next_node_cnt << ", "
         << next_edge_cnt << endl;
    if (max_n_adj - 1 <= level) break;
    level++;
    max_dim_s *= max_dim;
  }
  CUDA_CHECK(cudaMemcpy(sepsets.data(), sepsets_d, size_sepsets,
                        cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaFree(G_d));
  CUDA_CHECK(cudaFree(sepsets_d));
  // stage 2: orient edges
  PDAG G_pdag;
  G_pdag.g = vector<vector<bool>>(n_node, vector<bool>(n_node));
  for (int i = 0; i < n_node; i++) {
    for (int j = 0; j < n_node; j++) {
      G_pdag.g.at(i).at(j) = G[i * n_node + j];
    }
  }
  orientation(G_pdag, sepsets);
  return G_pdag;
}

py::array_t<bool> cuPC(py::array_t<uint8_t> data, py::array_t<int> n_states) {
  // translate input data to c++ vector(this is not optimal but I don't know
  // how to use pybind11::array_t)
  py::buffer_info buf_data = data.request(), buf_states = n_states.request();
  const uint8_t *__restrict__ prt_data = static_cast<uint8_t *>(buf_data.ptr);
  const int *__restrict__ prt_states = static_cast<int *>(buf_states.ptr);
  size_t n_data = buf_data.shape[0],
         n_node = buf_data.shape[1];  // number of nodes
  cout << "n_data, n_node: " << n_data << ' ' << n_node << endl;
  vector<uint8_t> data_vec(n_data * n_node);
  vector<int> n_states_vec(n_node);
  for (size_t i = 0; i < n_data; i++) {
    for (size_t j = 0; j < n_node; j++) {
      data_vec.at(j * n_data + i) = prt_data[i * n_node + j];
    }
  }
  for (size_t i = 0; i < n_node; i++) {
    n_states_vec.at(i) = prt_states[i];
  }
  for (size_t i = 0; i < n_data; i++) {
    for (size_t j = 0; j < n_node; j++) {
      assert(data_vec.at(j * n_data + i) < n_states_vec.at(j));
    }
  }
  auto endg = py::array_t<bool>({n_node, n_node});
  py::buffer_info buf_endg = endg.request();
  bool *__restrict__ prt_endg = static_cast<bool *>(buf_endg.ptr);

  PDAG Gend = PCsearch(n_node, n_data, data_vec, n_states_vec);

  // translate Gend to py::array_t (this is not optimal but I don't know how
  // to use pybind11::array_t)
  for (size_t i = 0; i < n_node; i++) {
    for (size_t j = 0; j < n_node; j++) {
      prt_endg[i * n_node + j] = Gend.g.at(i).at(j);
    }
  }
  return endg;
}
}  // namespace cuda3