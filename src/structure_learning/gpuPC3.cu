#include <cassert>
#include <cstddef>
#include <set>
#include <stack>
#include <string>
#include <vector>
using namespace std;

#include "structure_learning/gpuPC3.h"
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
// Based on: "GPU-Accelerated Constraint-Based Causal Structure Learning for
// Discrete Data.", 2021
//   authors: Hagedorn, Christopher, and Johannes Huegle
//   journal: Proceedings of the 2021 SIAM International Conference on Data
//   Mining (SDM)

// poz is a helper function for the afterwards defined pchisq
// It calculates the probability of a normal z score
// Implementation of the following function is adapted from Gary Perlman,
// source retrievable at
// https://www.netlib.org/a/perlman
//    Module:       z.c
//    Purpose:      compute approximations to normal z distribution
//    probabilities Programmer:   Gary Perlman Organization: Wang Institute,
//    Tyngsboro, MA 01879 Tester:       compile with -DZTEST to include main
//    program Copyright:    none Tabstops:     4
//
// z: z-score
#define Z_MAX 6.0                               /* maximum meaningful z value */
#define LOG_SQRT_PI 0.5723649429247000870717135 /* log (sqrt (pi)) */
#define I_SQRT_PI 0.5641895835477562869480795   /* 1 / sqrt (pi) */
#define BIGX 20.0 /* max value to represent exp (x) */
#define ex(x) (((x) < -BIGX) ? 0.0 : exp(x))
__device__ double poz(double z) {
  double y;
  double x;
  double w;

  if (z == 0.0) {
    x = 0.0;
  } else {
    y = 0.5 * fabs(z);
    if (y >= (Z_MAX * 0.5)) {
      x = 1.0;
    } else if (y < 1.0) {
      w = y * y;
      x = ((((((((0.000124818987 * w - 0.001075204047) * w + 0.005198775019) *
                    w -
                0.019198292004) *
                   w +
               0.059054035642) *
                  w -
              0.151968751364) *
                 w +
             0.319152932694) *
                w -
            0.531923007300) *
               w +
           0.797884560593) *
          y * 2.0;
    } else {
      y -= 2.0;
      x = (((((((((((((-0.000045255659 * y + 0.000152529290) * y -
                      0.000019538132) *
                         y -
                     0.000676904986) *
                        y +
                    0.001390604284) *
                       y -
                   0.000794620820) *
                      y -
                  0.002034254874) *
                     y +
                 0.006549791214) *
                    y -
                0.010557625006) *
                   y +
               0.011630447319) *
                  y -
              0.009279453341) *
                 y +
             0.005353579108) *
                y -
            0.002141268741) *
               y +
           0.000535310849) *
              y +
          0.999936657524;
    }
  }
  return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}

// Calculates the probability of a given chi-squared value
// x: obtained chi-squared value
// df: degrees of freedom

// Implementation of the following function
// is adapted from Gary Perlman, source retrievable at
// https://www.netlib.org/a/perlman
//    Module:       chisq.c
//    Purpose:      compute approximations to chisquare distribution
//    probabilities Contents:     pochisq(), critchi() Uses:         poz() in
//    z.c (Algorithm 209) Programmer:   Gary Perlman Organization: Wang
//    Institute, Tyngsboro, MA 01879 Tester:       compile with -DCHISQTEST to
//    include main program Copyright:    none Tabstops:     4
// which itself adapted from:
//            Hill, I. D. and Pike, M. C.  Algorithm 299
//            Collected Algorithms for the CACM 1967 p. 243
//    Updated for rounding errors based on remark in
//            ACM TOMS June 1985, page 185
__device__ double pchisq(double x, int df) {
  double a;
  double y = 0.0;
  double s;

  double e;
  double c;
  double z;

  bool even; /* true if df is an even number */

  if (x <= 0.0 || df < 1) {
    return (1.0);
  }

  a = 0.5 * x;
  even = (2 * (df / 2)) == df;
  if (df > 1) {
    y = ex(-a);
  }
  s = (even ? y : (2.0 * poz(-sqrt(x))));

  if (df <= 2) {
    return (s);
  }

  x = 0.5 * (df - 1.0);
  z = (even ? 1.0 : 0.5);
  if (a > BIGX) {
    e = (even ? 0.0 : LOG_SQRT_PI);
    c = log(a);
    while (z <= x) {
      e = log(z) + e;
      s += ex(c * z - a - e);
      z += 1.0;
    }
    return (s);
  }

  e = (even ? 1.0 : (I_SQRT_PI / sqrt(a)));
  c = 0.0;
  while (z <= x) {
    e = e * (a / z);
    c = c + e;
    z += 1.0;
  }
  return (c * y + s);
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

const int max_level = 20;
const int max_dim = 21;

__device__ bool ci_test_chi_squared_level_0(int n_data, int n_i, int n_j,
                                            int *contingency_matrix,
                                            int *marginals_i,
                                            int *marginals_j) {
  double chi_squared = 0;
  for (int k = 0; k < n_i; k++) {
    for (int l = 0; l < n_j; l++) {
      double expected =
          static_cast<double>(marginals_i[k]) * marginals_j[l] / n_data;
      if (expected != 0) {
        double observed = contingency_matrix[k * n_j + l];
        chi_squared += (observed - expected) * (observed - expected) / expected;
      }
    }
  }
  double pval = pchisq(chi_squared, (n_i - 1) * (n_j - 1));
  return pval >= 0.01;
}

__device__ bool ci_test_mi_level_0(int n_data, int n_i, int n_j,
                                   int *contingency_matrix, int *marginals_i,
                                   int *marginals_j) {
  double mi = 0;
  for (int k = 0; k < n_i; k++) {
    for (int l = 0; l < n_j; l++) {
      if (!contingency_matrix[k * n_j + l]) continue;
      mi += static_cast<double>(contingency_matrix[k * n_j + l]) / n_data *
            log2(static_cast<double>(n_data) * contingency_matrix[k * n_j + l] /
                 (static_cast<double>(marginals_i[k]) * marginals_j[l]));
    }
  }
  return mi < 0.003;
}

__device__ bool ci_test_bayes_factor_level_0(int n_data, int n_i, int n_j,
                                             int *contingency_matrix,
                                             int *marginals_i,
                                             int *marginals_j) {
  double independent_score = 0;
  double dependent_score = 0;
  const double alpha = 0.5;
  for (int k = 0; k < n_i; k++) {
    independent_score += lgamma(marginals_i[k] + alpha) - lgamma(alpha);
  }
  independent_score += lgamma(n_i * alpha) - lgamma(n_i * alpha + n_data);
  for (int l = 0; l < n_j; l++) {
    independent_score += lgamma(marginals_j[l] + alpha) - lgamma(alpha);
  }
  independent_score += lgamma(n_j * alpha) - lgamma(n_j * alpha + n_data);
  for (int k = 0; k < n_i; k++) {
    for (int l = 0; l < n_j; l++) {
      dependent_score +=
          lgamma(contingency_matrix[k * n_j + l] + alpha) - lgamma(alpha);
    }
  }
  dependent_score +=
      lgamma(n_i * n_j * alpha) - lgamma(n_i * n_j * alpha + n_data);
  return independent_score > dependent_score - 1e-10;
}

__global__ void PC_level_0(int n_node, int n_data, uint8_t *data, int *G,
                           int *n_states) {
  int i = blockIdx.x;
  int j = blockIdx.y;
  if (i >= j) {
    return;
  }
  __shared__ int contingency_matrix[max_dim * max_dim];
  int n_i = n_states[i];
  int n_j = n_states[j];
  for (int k = threadIdx.x; k < n_i * n_j; k += blockDim.x) {
    contingency_matrix[k] = 0;
  }
  __syncthreads();
  for (int k = threadIdx.x; k < n_data; k += blockDim.x) {
    int idx = data[i * n_data + k] * n_j + data[j * n_data + k];
    atomicAdd(contingency_matrix + idx, 1);
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    int marginals_i[max_dim];
    int marginals_j[max_dim];
    for (int k = 0; k < max_dim; k++) {
      marginals_i[k] = 0;
      marginals_j[k] = 0;
    }
    for (int k = 0; k < n_i; k++) {
      for (int l = 0; l < n_j; l++) {
        int entry = contingency_matrix[k * n_j + l];
        marginals_i[k] += entry;
        marginals_j[l] += entry;
      }
    }
    if (ci_test_chi_squared_level_0(n_data, n_i, n_j, contingency_matrix,
                                    marginals_i, marginals_j)) {
      G[i * n_node + j] = 0;
      G[j * n_node + i] = 0;
    }
  }
}

__device__ void ci_test_chi_squared_level_n(double *chi_squared, int n_data,
                                            int dim_s, int dim_mul, int n_i,
                                            int n_j, int *N_i_j_s, int *N_i_s,
                                            int *N_j_s, int *N_s,
                                            bool *result) {
  if (threadIdx.x == 0) {
    *chi_squared = 0;
  }
  __syncthreads();
  for (int g = threadIdx.x; g < dim_s; g += blockDim.x) {
    int h = g / dim_mul / n_j * dim_mul + g % dim_mul;
    if (N_s[h] == 0) continue;
    int l = (g / dim_mul) % n_j;
    for (int k = 0; k < n_i; k++) {
      double expected =
          static_cast<double>(N_i_s[h * n_i + k]) * N_j_s[h * n_j + l] / N_s[h];
      if (expected == 0) continue;
      double observed = N_i_j_s[g * n_i + k];
      double sum_term =
          (observed - expected) * (observed - expected) / expected;
      unsigned long long *address_as_ull =
          reinterpret_cast<unsigned long long *>(chi_squared);
      unsigned long long old = *address_as_ull, assumed;
      do {
        assumed = old;
        old = atomicCAS(
            address_as_ull, assumed,
            __double_as_longlong(sum_term + __longlong_as_double(assumed)));
      } while (assumed != old);
    }
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    double pval = pchisq(*chi_squared, (n_i - 1) * (n_j - 1) * dim_s / n_j);
    *result = (pval >= 0.01);
  }
}

__device__ void ci_test_mi_level_n(double *mi, int n_data, int dim_s, int n_i,
                                   int n_j, int *N_i_j_s, int *N_i_s,
                                   int *N_j_s, int *N_s, bool *result) {
  if (threadIdx.x == 0) {
    *mi = 0;
  }
  __syncthreads();
  for (int g = threadIdx.x; g < dim_s; g += blockDim.x) {
    if (N_s[g] == 0) continue;
    for (int k = 0; k < n_i; k++) {
      for (int l = 0; l < n_j; l++) {
        if (!N_i_j_s[g * n_i * n_j + k * n_j + l]) continue;
        double sum_term =
            static_cast<double>(N_i_j_s[g * n_i * n_j + k * n_j + l]) / n_data *
            log2(
                static_cast<double>(N_s[g]) *
                N_i_j_s[g * n_i * n_j + k * n_j + l] /
                (static_cast<double>(N_i_s[g * n_i + k]) * N_j_s[g * n_j + l]));
        unsigned long long *address_as_ull =
            reinterpret_cast<unsigned long long *>(mi);
        unsigned long long old = *address_as_ull, assumed;
        do {
          assumed = old;
          old = atomicCAS(
              address_as_ull, assumed,
              __double_as_longlong(sum_term + __longlong_as_double(assumed)));
        } while (assumed != old);
      }
    }
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    *result = (*mi < 0.003);
  }
}

__device__ void ci_test_bayes_factor_level_n(double *scratch_ptr, int n_data,
                                             int dim_s, int n_i, int n_j,
                                             int *N_i_j_s, int *N_i_s,
                                             int *N_j_s, int *N_s,
                                             bool *result) {
  if (threadIdx.x == 0) {
    *scratch_ptr = 0;
  }
  __syncthreads();
  const double alpha = 0.5;
  // independent score
  for (int g = threadIdx.x; g < dim_s; g += blockDim.x) {
    if (N_s[g] == 0) continue;
    double sum_term = 0;
    for (int k = 0; k < n_i; k++) {
      sum_term += lgamma(N_i_s[g * n_i + k] + alpha) - lgamma(alpha);
    }
    sum_term += lgamma(n_i * alpha) - lgamma(n_i * alpha + N_s[g]);
    for (int l = 0; l < n_j; l++) {
      sum_term += lgamma(N_j_s[g * n_j + l] + alpha) - lgamma(alpha);
    }
    sum_term += lgamma(n_j * alpha) - lgamma(n_j * alpha + N_s[g]);
    unsigned long long *address_as_ull =
        reinterpret_cast<unsigned long long *>(scratch_ptr);
    unsigned long long old = *address_as_ull, assumed;
    do {
      assumed = old;
      old = atomicCAS(
          address_as_ull, assumed,
          __double_as_longlong(sum_term + __longlong_as_double(assumed)));
    } while (assumed != old);
  }
  __syncthreads();
  double independent_score = *scratch_ptr;
  // dependent score
  if (threadIdx.x == 0) {
    *scratch_ptr = 0;
  }
  __syncthreads();
  for (int g = threadIdx.x; g < dim_s; g += blockDim.x) {
    if (N_s[g] == 0) continue;
    double sum_term = 0;
    for (int k = 0; k < n_i; k++) {
      for (int l = 0; l < n_j; l++) {
        sum_term += lgamma(N_i_j_s[g * n_i * n_j + k * n_j + l] + alpha) -
                    lgamma(alpha);
      }
    }
    sum_term += lgamma(n_i * n_j * alpha) - lgamma(n_i * n_j * alpha + N_s[g]);
    unsigned long long *address_as_ull =
        reinterpret_cast<unsigned long long *>(scratch_ptr);
    unsigned long long old = *address_as_ull, assumed;
    do {
      assumed = old;
      old = atomicCAS(
          address_as_ull, assumed,
          __double_as_longlong(sum_term + __longlong_as_double(assumed)));
    } while (assumed != old);
  }
  __syncthreads();
  double dependent_score = *scratch_ptr;
  if (threadIdx.x == 0) {
    // printf("independent_score: %.7lf, dependent_score: %.7lf\n",
    //        independent_score, dependent_score);
    *result = (independent_score > dependent_score - 1e-10);
  }
}

__global__ void PC_level_n(int level, int i_offset, int n_node, int n_data,
                           uint8_t *data, int *G, int *n_states,
                           int *working_memory, int *sepsets) {
  extern __shared__ int smem[];
  int i = i_offset + blockIdx.x;
  if (i >= n_node) return;
  int *G_compacted = smem;
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    int cnt = 0;
    for (int j = 0; j < n_node; j++) {
      if (G[i * n_node + j]) {
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
  int max_dim_s = pow(static_cast<double>(max_dim), level);
  int reserved_size_per_ci_test =
      max_dim_s * max_dim * max_dim + 2 * max_dim_s * max_dim + max_dim_s;
  int sepset_cnt = binom(n_adj, level + 1);
  int step = gridDim.y * blockDim.y;
  int sepset[max_level + 1];
  int dim_mul[max_level + 1];
  int ci_test_idx = blockIdx.y * blockDim.y + threadIdx.y;
  int thread_memory_index = step * blockIdx.x + ci_test_idx;
  int *thread_memory =
      working_memory + thread_memory_index * reserved_size_per_ci_test;
  for (int sepset_idx = ci_test_idx; sepset_idx < sepset_cnt;
       sepset_idx += step) {
    comb(n_adj, level + 1, sepset_idx, -1, sepset);
    int dim_s = 1;
    bool valid = false;
    for (int k = 0; k < level + 1; k++) {
      sepset[k] = G_compacted[sepset[k] + 1];
      if (G[i * n_node + sepset[k]] == 1) valid = true;
      dim_s *= n_states[sepset[k]];
    }
    if (!valid) continue;
    __syncthreads();
    dim_mul[level] = 1;
    for (int k = level - 1; k >= 0; k--) {
      dim_mul[k] = dim_mul[k + 1] * n_states[sepset[k + 1]];
    }
    int *N_i_j_s = thread_memory;
    int n_i = n_states[i];
    if (threadIdx.x == 0) {
      int malloc_size = dim_s * n_i;
      memset(N_i_j_s, 0, malloc_size * sizeof(int));
    }
    __syncthreads();
    for (int k = threadIdx.x; k < n_data; k += blockDim.x) {
      int val_i = data[i * n_data + k];
      int s_idx = 0;
      for (int l = 0; l < level + 1; l++) {
        s_idx = s_idx * n_states[sepset[l]] + data[sepset[l] * n_data + k];
      }
      atomicAdd(N_i_j_s + s_idx * n_i + val_i, 1);
    }
    for (int idx_j = 0; idx_j < level + 1; idx_j++) {
      int j = sepset[idx_j];
      if (G[i * n_node + j] != 1) {
        continue;
      }
      int n_j = n_states[j];
      int dim_mul_j = dim_mul[idx_j];
      int *N_i_s = N_i_j_s + dim_s * n_i;
      int *N_j_s = N_i_s + dim_s * n_i / n_j;
      int *N_s = N_j_s + dim_s;
      if (threadIdx.x == 0) {
        int malloc_size = dim_s * n_i / n_j + dim_s + dim_s / n_j;
        memset(N_i_s, 0, malloc_size * sizeof(int));
      }
      __syncthreads();
      for (int g = threadIdx.x; g < dim_s; g += blockDim.x) {
        for (int k = 0; k < n_i; k++) {
          int l = (g / dim_mul_j) % n_j;
          int h = g / dim_mul_j / n_j * dim_mul_j + g % dim_mul_j;
          int entry = N_i_j_s[g * n_i + k];
          atomicAdd(N_i_s + h * n_i + k, entry);
          atomicAdd(N_j_s + h * n_j + l, entry);
          atomicAdd(N_s + h, entry);
        }
      }
      __syncthreads();
      int scratch_addr = n_adj + 1 + threadIdx.y;
      scratch_addr = (scratch_addr + 1) / 2 * 2;
      double *scratch_ptr = reinterpret_cast<double *>(smem + scratch_addr);
      bool result;
      ci_test_chi_squared_level_n(scratch_ptr, n_data, dim_s, dim_mul_j, n_i,
                                  n_j, N_i_j_s, N_i_s, N_j_s, N_s, &result);
      if (threadIdx.x == 0 && result) {
        if (atomicCAS(G + i * n_node + j, 1, -1) == 1) {
          G[j * n_node + i] = -1;
          int ij_min = (i < j ? i : j);
          int ij_max = (i < j ? j : i);
          int p = 0;
          for (int k = 0; k < level + 1; k++) {
            if (k == idx_j) continue;
            sepsets[(ij_min * n_node + ij_max) * max_level + p] = sepset[k];
            p++;
          }
        }
      }
      __syncthreads();
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
  vector<int> sepsets(n_node * n_node * max_level, -1);
  uint8_t *data_d;
  int *G_d, *n_states_d, *working_memory_d, *sepsets_d;
  int size_G = sizeof(int) * n_node * n_node;
  int size_data = sizeof(uint8_t) * n_data * n_node;
  int size_n_states = sizeof(int) * n_node;
  int size_working_memory = sizeof(int) * 500'000'000;
  int size_sepsets = sizeof(int) * n_node * n_node * max_level;
  CUDA_CHECK(cudaMalloc(&G_d, size_G));
  CUDA_CHECK(cudaMalloc(&data_d, size_data));
  CUDA_CHECK(cudaMalloc(&n_states_d, size_n_states));
  CUDA_CHECK(cudaMalloc(&working_memory_d, size_working_memory));
  CUDA_CHECK(cudaMalloc(&sepsets_d, size_sepsets));
  CUDA_CHECK(
      cudaMemcpy(data_d, data.data(), size_data, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(n_states_d, n_states.data(), size_n_states,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(sepsets_d, sepsets.data(), size_sepsets,
                        cudaMemcpyHostToDevice));

  // stage 1: Do CI tests between nodes and remove edge (undirected graph)
  int level = 0;
  int max_n_adj = n_node - 1;
  uint64_t max_dim_s = 1;
  while (level <= n_node - 2) {
    CUDA_CHECK(cudaMemcpy(G_d, G.data(), size_G, cudaMemcpyHostToDevice));
    cout << "level: " << level << ", max_n_adj: " << max_n_adj << endl;
    if (level == 0) {
      dim3 threadsPerBlock(64);
      dim3 numBlocks(n_node, n_node);
      PC_level_0<<<numBlocks, threadsPerBlock>>>(n_node, n_data, data_d, G_d,
                                                 n_states_d);
    } else {
      dim3 threadsPerBlock(64, 1);
      dim3 numBlocks(n_node, max_n_adj * 2);
      uint64_t reserved_size_per_ci_test =
          max_dim_s * max_dim * max_dim + 2 * max_dim_s * max_dim + max_dim_s;
      if (reserved_size_per_ci_test > size_working_memory / sizeof(int)) {
        cout << "working memory is not enough" << endl;
        break;
      }
      uint64_t reserved_size_per_row =
          reserved_size_per_ci_test * threadsPerBlock.y * numBlocks.y;
      int max_rows = size_working_memory / sizeof(int) / reserved_size_per_row;
      if (max_rows == 0) {
        int max_columns =
            size_working_memory / sizeof(int) / reserved_size_per_ci_test;
        numBlocks.x = 1;
        numBlocks.y = max_columns;
        threadsPerBlock.y = 1;
      } else if (numBlocks.x > max_rows) {
        numBlocks.x = max_rows;
      }
      cout << "numBlocks: " << numBlocks.x << ", " << numBlocks.y << endl;
      cout << "threadsPerBlock: " << threadsPerBlock.x << ", "
           << threadsPerBlock.y << endl;
      int loop = (n_node + numBlocks.x - 1) / numBlocks.x;
      for (int l = 0; l < loop; l++) {
        int i_offset = l * numBlocks.x;
        PC_level_n<<<numBlocks, threadsPerBlock,
                     sizeof(int) * (max_n_adj) +
                         sizeof(double) * (threadsPerBlock.y + 1)>>>(
            level, i_offset, n_node, n_data, data_d, G_d, n_states_d,
            working_memory_d, sepsets_d);
        // cudaDeviceSynchronize();
      }
    }
    CUDA_CHECK(cudaMemcpy(G.data(), G_d, size_G, cudaMemcpyDeviceToHost));
    max_n_adj = 0;
    int next_node_cnt = 0;
    int next_edge_cnt = 0;
    for (int i = 0; i < n_node; i++) {
      int n_adj = 0;
      for (int j = 0; j < n_node; j++) {
        if (G[i * n_node + j] == -1) {
          G[i * n_node + j] = 0;
        }
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
  CUDA_CHECK(cudaFree(data_d));
  CUDA_CHECK(cudaFree(n_states_d));
  CUDA_CHECK(cudaFree(working_memory_d));
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

py::array_t<bool> gpuPC3(py::array_t<uint8_t> data, py::array_t<int> n_states) {
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