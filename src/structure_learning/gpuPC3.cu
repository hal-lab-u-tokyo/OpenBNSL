#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace std;

#include "base/PDAG2.h"
#include "structure_learning/constants.h"
#include "structure_learning/gpuPC3.h"
#include "structure_learning/orientation.h"

namespace cuda3 {
#include "structure_learning/utils.cuh"

__global__ void calc_regret(double *regret) {
  int n = blockDim.x * blockIdx.x + threadIdx.x;
  if (n == 0) return;
  double sum = 1, a = 1;
  int bound = ceil(2 + sqrt(20 * n * log(10.0)));
  for (int k = 1; k <= bound; k++) {
    a *= (n - k + 1) / static_cast<double>(n);
    sum += a;
  }
  regret[n * max_dim] = 1;
  regret[n * max_dim + 1] = sum;
  for (int k = 3; k <= max_dim; k++) {
    regret[n * max_dim + k - 1] =
        regret[n * max_dim + k - 2] +
        n * regret[n * max_dim + k - 3] / static_cast<double>(k - 2);
  }
  for (int k = 0; k < max_dim; k++) {
    regret[n * max_dim + k] = log2(max(1.0, regret[n * max_dim + k]));
  }
}

__device__ bool ci_test_sc_level_0(int n_data, int n_i, int n_j,
                                   int *contingency_matrix, int *marginals_i,
                                   int *marginals_j, double *regret) {
  double mi = 0;
  for (int k = 0; k < n_i; k++) {
    for (int l = 0; l < n_j; l++) {
      if (!contingency_matrix[k * n_j + l]) continue;
      mi += static_cast<double>(contingency_matrix[k * n_j + l]) / n_data *
            log2(static_cast<double>(n_data) * contingency_matrix[k * n_j + l] /
                 (static_cast<double>(marginals_i[k]) * marginals_j[l]));
    }
  }
  // R(X_i)
  double r_i = regret[n_data * max_dim + n_i - 1];
  // R(X_i|X_j)
  double r_ij = 0;
  for (int l = 0; l < n_j; l++) {
    r_ij += regret[marginals_j[l] * max_dim + n_i - 1];
  }
  // R(X_j)
  double r_j = regret[n_data * max_dim + n_j - 1];
  double r_ji = 0;
  for (int l = 0; l < n_i; l++) {
    r_ji += regret[marginals_i[l] * max_dim + n_j - 1];
  }
  double threshold = min(r_ij - r_i, r_ji - r_j) / n_data;
  // printf("threshold: %.7lf\n", threshold);
  return mi <= threshold;
}

__device__ bool ci_test_g2_level_0(int n_data, int n_i, int n_j,
                                   int *contingency_matrix, int *marginals_i,
                                   int *marginals_j) {
  if (n_i == 1 || n_j == 1) {
    return true;
  }
  double g2 = 0;
  for (int k = 0; k < n_i; k++) {
    for (int l = 0; l < n_j; l++) {
      double expected =
          static_cast<double>(marginals_i[k]) * marginals_j[l] / n_data;
      double observed = contingency_matrix[k * n_j + l];
      if (observed != 0) {
        g2 += 2 * observed * log(observed / expected);
      }
    }
  }
  double pval = pchisq(g2, (n_i - 1) * (n_j - 1));
  return pval >= 0.05;
}

__global__ void PC_level_0(int citest_type, int n_node, int n_data,
                           uint8_t *data, int *G, int *n_states, double *regret,
                           int *model, int *stats) {
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
  if (threadIdx.x == 0) {
    uint smid;
    asm volatile("mov.u32 %0, %smid;" : "=r"(smid));
    atomicAdd(stats + smid, 1);
    atomicAdd(stats + smid + sm_num, 1);
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
    bool result;
    if (citest_type == 0) {
      result = ci_test_g2_level_0(n_data, n_i, n_j, contingency_matrix,
                                  marginals_i, marginals_j);
    } else {
      result = ci_test_sc_level_0(n_data, n_i, n_j, contingency_matrix,
                                  marginals_i, marginals_j, regret);
    }
    if (result) {
      G[i * n_node + j] = 0;
      G[j * n_node + i] = 0;
    }
    // bool sep_result = d_separated(0, n_node, i, j, nullptr, model);
    // if (result && sep_result) {
    //   atomicAdd(stats, 1);
    // } else if (result && !sep_result) {
    //   atomicAdd(stats + 1, 1);
    // } else if (!result && sep_result) {
    //   atomicAdd(stats + 2, 1);
    // } else {
    //   atomicAdd(stats + 3, 1);
    // }
  }
}

__device__ void ci_test_sc_level_n(double *mi, int n_data, int dim_s,
                                   int dim_mul, int n_i, int n_j, int *N_i_j_s,
                                   int *N_i_s, int *N_j_s, int *N_s,
                                   bool *result, double *regret) {
  double *r_i = mi + 1;
  double *r_ij = mi + 2;
  double *r_j = mi + 3;
  double *r_ji = mi + 4;
  if (threadIdx.x == 0) {
    *mi = *r_i = *r_ij = *r_j = *r_ji = 0;
  }
  __syncthreads();
  for (int h = threadIdx.x; h < dim_s; h += blockDim.x) {
    if (N_s[h] == 0) continue;
    atomicAdd(r_i, regret[N_s[h] * max_dim + n_i - 1]);
    atomicAdd(r_j, regret[N_s[h] * max_dim + n_j - 1]);
    for (int k = 0; k < n_i; k++) {
      atomicAdd(r_ji, regret[N_i_s[h * n_i + k] * max_dim + n_j - 1]);
      if (k && !N_i_s[h * n_i + k]) continue;
      for (int l = 0; l < n_j; l++) {
        if (k == 0) {
          atomicAdd(r_ij, regret[N_j_s[h * n_j + l] * max_dim + n_i - 1]);
        }
        int g = (h / dim_mul) * dim_mul * n_j + l * dim_mul + h % dim_mul;
        double entry = N_i_j_s[g * n_i + k];
        if (entry) {
          atomicAdd(mi, entry / n_data *
                            log2(N_s[h] * entry /
                                 (static_cast<double>(N_i_s[h * n_i + k]) *
                                  N_j_s[h * n_j + l])));
        }
      }
    }
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    double threshold = min(*r_ij - *r_i, *r_ji - *r_j) / n_data;
    *result = (*mi <= threshold);
  }
}

__device__ void ci_test_sc_level_n_2(double *mi, int n_data, int dim_s,
                                     int dim_mul_i, int dim_mul_j, int n_i,
                                     int n_j, int *N_i_j_s, int *N_i_s,
                                     int *N_j_s, int *N_s, bool *result,
                                     double *regret) {
  double *r_i = mi + 1;
  double *r_ij = mi + 2;
  double *r_j = mi + 3;
  double *r_ji = mi + 4;
  if (threadIdx.x == 0) {
    *mi = *r_i = *r_ij = *r_j = *r_ji = 0;
  }
  __syncthreads();
  for (int h = threadIdx.x; h < dim_s; h += blockDim.x) {
    if (N_s[h] == 0) continue;
    atomicAdd(r_i, regret[N_s[h] * max_dim + n_i - 1]);
    atomicAdd(r_j, regret[N_s[h] * max_dim + n_j - 1]);
    for (int k = 0; k < n_i; k++) {
      atomicAdd(r_ji, regret[N_i_s[h * n_i + k] * max_dim + n_j - 1]);
      for (int l = 0; l < n_j; l++) {
        if (k == 0) {
          atomicAdd(r_ij, regret[N_j_s[h * n_j + l] * max_dim + n_i - 1]);
        }
        int g = (h / (dim_mul_j / n_i)) * dim_mul_j * n_j + l * dim_mul_j +
                (h % (dim_mul_j / n_i)) / dim_mul_i * dim_mul_i * n_i +
                k * dim_mul_i + h % dim_mul_i;
        double entry = N_i_j_s[g];
        if (entry) {
          atomicAdd(mi, entry / n_data *
                            log2(N_s[h] * entry /
                                 (static_cast<double>(N_i_s[h * n_i + k]) *
                                  N_j_s[h * n_j + l])));
        }
      }
    }
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    double threshold = min(*r_ij - *r_i, *r_ji - *r_j) / n_data;
    *result = (*mi <= threshold);
  }
}

__device__ void ci_test_g2_level_n(double *g2, int n_data, int dim_s,
                                   int dim_mul, int n_i, int n_j, int *N_i_j_s,
                                   int *N_i_s, int *N_j_s, int *N_s,
                                   bool *result) {
  int *df = reinterpret_cast<int *>(g2 + 1);
  if (threadIdx.x == 0) {
    *g2 = 0;
    *df = 0;
  }
  __syncthreads();
  for (int h = threadIdx.x; h < dim_s; h += blockDim.x) {
    if (N_s[h] == 0) continue;
    int alx = 0, aly = 0;
    for (int l = 0; l < n_j; l++) {
      aly += (N_j_s[h * n_j + l] > 0);
    }
    for (int k = 0; k < n_i; k++) {
      if (!N_i_s[h * n_i + k]) continue;
      alx++;
      for (int l = 0; l < n_j; l++) {
        int g = (h / dim_mul) * dim_mul * n_j + l * dim_mul + h % dim_mul;
        double expected = static_cast<double>(N_i_s[h * n_i + k]) *
                          N_j_s[h * n_j + l] / N_s[h];
        double observed = N_i_j_s[g * n_i + k];
        if (observed) {
          double sum_term = 2 * observed * log(observed / expected);
          atomicAdd(g2, sum_term);
        }
      }
    }
    alx = (alx >= 1 ? alx : 1);
    aly = (aly >= 1 ? aly : 1);
    atomicAdd(df, (alx - 1) * (aly - 1));
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    if (df == 0) {
      *result = true;
    } else {
      double pval = pchisq(*g2, *df);
      *result = (pval >= 0.05);
    }
  }
}

__device__ void ci_test_g2_level_n_2(double *g2, int n_data, int dim_s,
                                     int dim_mul_i, int dim_mul_j, int n_i,
                                     int n_j, int *N_i_j_s, int *N_i_s,
                                     int *N_j_s, int *N_s, bool *result) {
  int *df = reinterpret_cast<int *>(g2 + 1);
  if (threadIdx.x == 0) {
    *g2 = 0;
    *df = 0;
  }
  __syncthreads();
  for (int h = threadIdx.x; h < dim_s; h += blockDim.x) {
    if (N_s[h] == 0) continue;
    int alx = 0, aly = 0;
    for (int l = 0; l < n_j; l++) {
      aly += (N_j_s[h * n_j + l] > 0);
    }
    for (int k = 0; k < n_i; k++) {
      if (!N_i_s[h * n_i + k]) continue;
      alx++;
      for (int l = 0; l < n_j; l++) {
        int g = (h / (dim_mul_j / n_i)) * dim_mul_j * n_j + l * dim_mul_j +
                (h % (dim_mul_j / n_i)) / dim_mul_i * dim_mul_i * n_i +
                k * dim_mul_i + h % dim_mul_i;
        double expected = static_cast<double>(N_i_s[h * n_i + k]) *
                          N_j_s[h * n_j + l] / N_s[h];
        double observed = N_i_j_s[g];
        if (observed) {
          double sum_term = 2 * observed * log(observed / expected);
          atomicAdd(g2, sum_term);
        }
      }
    }
    alx = (alx >= 1 ? alx : 1);
    aly = (aly >= 1 ? aly : 1);
    atomicAdd(df, (alx - 1) * (aly - 1));
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    if (df == 0) {
      *result = true;
    } else {
      double pval = pchisq(*g2, *df);
      *result = (pval >= 0.05);
    }
  }
}

__global__ void PC_level_n(int citest_type, int level, int n_node, int n_data,
                           uint8_t *data, int *G, int *n_states,
                           bool use_working_memory, int *working_memory,
                           int *sepsets, double *regret, int *model,
                           int *stats) {
  extern __shared__ int smem[];
  for (int i = blockIdx.x; i < n_node; i += gridDim.x) {
    __syncthreads();
    int *G_compacted = smem;
    if (threadIdx.x == 0) {
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
      continue;
    }
    int max_dim_s = pow(static_cast<double>(max_dim), level);
    int reserved_size_per_ci_test =
        max_dim_s * max_dim * max_dim + 2 * max_dim_s * max_dim + max_dim_s;
    int sepset_cnt = binom(n_adj, level + 1);
    int *sepset = smem + n_adj + 2;
    int *dim_mul = smem + n_adj + 2 + level + 1;
    int *connected_to_all = smem + n_adj + 2 + (level + 1) * 2;
    int *thread_memory;
    if (use_working_memory) {
      thread_memory = working_memory + (gridDim.y * blockIdx.x + blockIdx.y) *
                                           reserved_size_per_ci_test;
    } else {
      thread_memory = smem + n_adj + 2 + (level + 1) * 3;
    }
    int n_i = n_states[i];
    for (int sepset_idx = blockIdx.y; sepset_idx < sepset_cnt;
         sepset_idx += gridDim.y) {
      __syncthreads();
      int *valid = smem + n_adj + 1;
      if (threadIdx.x == 0) {
        comb(n_adj, level + 1, sepset_idx, -1, sepset);
        *valid = 0;
        dim_mul[0] = 1;
        for (int k = 0; k < level + 1; k++) {
          sepset[k] = G_compacted[sepset[k] + 1];
          if (G[i * n_node + sepset[k]] == 1) {
            *valid = 1;
          }
          if (k < level) {
            dim_mul[k + 1] = dim_mul[k] * n_states[sepset[k]];
          }
        }
        for (int idx_j = 0; idx_j < level + 1; idx_j++) {
          int j = sepset[idx_j];
          connected_to_all[idx_j] = 1;
          bool exist_nondeleted_edge = false;
          for (int idx_k = 0; idx_k < level + 1; idx_k++) {
            if (idx_j != idx_k && !G[j * n_node + sepset[idx_k]]) {
              connected_to_all[idx_j] = 0;
              break;
            }
            if (idx_j != idx_k && G[j * n_node + sepset[idx_k]] == 1) {
              exist_nondeleted_edge = true;
            }
          }
          if (connected_to_all[idx_j]) {
            if (j < i) {
              *valid = 0;
              break;
            }
            if (j > i && exist_nondeleted_edge) {
              *valid = 1;
            }
          }
        }
      }
      __syncthreads();
      if (*valid == 0) continue;
      if (threadIdx.x == 0) {
        uint smid;
        asm volatile("mov.u32 %0, %smid;" : "=r"(smid));
        atomicAdd(stats + smid + sm_num, 1);
      }
      int dim_s = dim_mul[level] * n_states[sepset[level]];
      int *N_i_j_s = thread_memory;
      if (threadIdx.x == 0) {
        int malloc_size = dim_s * n_i;
        memset(N_i_j_s, 0, malloc_size * sizeof(int));
      }
      __syncthreads();
      for (int k = threadIdx.x; k < n_data; k += blockDim.x) {
        int val_i = data[i * n_data + k];
        int s_idx = 0;
        for (int l = 0; l < level + 1; l++) {
          s_idx += data[sepset[l] * n_data + k] * dim_mul[l];
        }
        atomicAdd(N_i_j_s + s_idx * n_i + val_i, 1);
      }
      int scratch_addr = n_adj + 2 + (level + 1) * 3 +
                         (use_working_memory ? 0 : reserved_size_per_ci_test);
      scratch_addr = (scratch_addr + 1) / 2 * 2;
      double *scratch_ptr = reinterpret_cast<double *>(smem + scratch_addr);
      bool result;
      for (int idx_j = 0; idx_j < level; idx_j++) {
        int j = sepset[idx_j];
        for (int idx_k = idx_j + 1; idx_k < level + 1; idx_k++) {
          int k = sepset[idx_k];
          if (threadIdx.x == 0) {
            if (G[j * n_node + k] != 1) {
              *valid = 0;
            } else {
              *valid = (connected_to_all[idx_j] || connected_to_all[idx_k]);
            }
          }
          __syncthreads();
          if (*valid == 0) continue;
          if (threadIdx.x == 0) {
            uint smid;
            asm volatile("mov.u32 %0, %smid;" : "=r"(smid));
            atomicAdd(stats + smid, 1);
          }
          int n_j = n_states[j];
          int n_k = n_states[k];
          int *N_i_s = N_i_j_s + dim_s * n_i;
          int *N_j_s = N_i_s + dim_s * n_i / n_k;
          int *N_s = N_j_s + dim_s * n_i / n_j;
          if (threadIdx.x == 0) {
            int malloc_size =
                dim_s * n_i / n_k + dim_s * n_i / n_j + dim_s * n_i / n_j / n_k;
            memset(N_i_s, 0, malloc_size * sizeof(int));
          }
          __syncthreads();
          for (int g = threadIdx.x; g < dim_s * n_i; g += blockDim.x) {
            int h =
                g / (dim_mul[idx_k] * n_i) / n_k * dim_mul[idx_k] * n_i / n_j +
                g % (dim_mul[idx_k] * n_i) / (dim_mul[idx_j] * n_i) / n_j *
                    dim_mul[idx_j] * n_i +
                g % (dim_mul[idx_j] * n_i);
            int k = g / (dim_mul[idx_j] * n_i) % n_j;
            int l = g / (dim_mul[idx_k] * n_i) % n_k;
            int entry = N_i_j_s[g];
            atomicAdd(N_i_s + h * n_j + k, entry);
            atomicAdd(N_j_s + h * n_k + l, entry);
            atomicAdd(N_s + h, entry);
          }
          if (citest_type == 0) {
            ci_test_g2_level_n_2(scratch_ptr, n_data, dim_s / n_j / n_k * n_i,
                                 dim_mul[idx_j] * n_i, dim_mul[idx_k] * n_i,
                                 n_j, n_k, N_i_j_s, N_i_s, N_j_s, N_s, &result);
          } else {
            ci_test_sc_level_n_2(scratch_ptr, n_data, dim_s / n_j / n_k * n_i,
                                 dim_mul[idx_j] * n_i, dim_mul[idx_k] * n_i,
                                 n_j, n_k, N_i_j_s, N_i_s, N_j_s, N_s, &result,
                                 regret);
          }
          if (threadIdx.x == 0 && result) {
            if (atomicCAS(G + j * n_node + k, 1, -1) == 1) {
              G[k * n_node + j] = -1;
              sepsets[(j * n_node + k) * max_level] = i;
              int p = 1;
              for (int l = 0; l < level + 1; l++) {
                if (l == idx_j || l == idx_k) continue;
                sepsets[(j * n_node + k) * max_level + p] = sepset[l];
                p++;
              }
            }
          }
        }
      }
      for (int idx_j = 0; idx_j < level + 1; idx_j++) {
        __syncthreads();
        int j = sepset[idx_j];
        if (threadIdx.x == 0) {
          *valid = (G[i * n_node + j] == 1);
        }
        __syncthreads();
        if (*valid == 0) continue;
        if (threadIdx.x == 0) {
          uint smid;
          asm volatile("mov.u32 %0, %smid;" : "=r"(smid));
          atomicAdd(stats + smid, 1);
        }
        int n_j = n_states[j];
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
            int h = (g / dim_mul[idx_j]) / n_j * dim_mul[idx_j] +
                    g % dim_mul[idx_j];
            int l = g / dim_mul[idx_j] % n_j;
            int entry = N_i_j_s[g * n_i + k];
            atomicAdd(N_i_s + h * n_i + k, entry);
            atomicAdd(N_j_s + h * n_j + l, entry);
            atomicAdd(N_s + h, entry);
          }
        }
        if (citest_type == 0) {
          ci_test_g2_level_n(scratch_ptr, n_data, dim_s / n_j, dim_mul[idx_j],
                             n_i, n_j, N_i_j_s, N_i_s, N_j_s, N_s, &result);
        } else {
          ci_test_sc_level_n(scratch_ptr, n_data, dim_s / n_j, dim_mul[idx_j],
                             n_i, n_j, N_i_j_s, N_i_s, N_j_s, N_s, &result,
                             regret);
        }
        if (threadIdx.x == 0 && result) {
          int ij_min = (i < j ? i : j);
          int ij_max = (i < j ? j : i);
          if (atomicCAS(G + ij_min * n_node + ij_max, 1, -1) == 1) {
            G[ij_max * n_node + ij_min] = -1;
            int p = 0;
            for (int k = 0; k < level + 1; k++) {
              if (k == idx_j) continue;
              sepsets[(ij_min * n_node + ij_max) * max_level + p] = sepset[k];
              p++;
            }
          }
        }
        // if (threadIdx.x == 0) {
        //   int sepset2[max_level];
        //   int p = 0;
        //   for (int k = 0; k < level + 1; k++) {
        //     if (k == idx_j) continue;
        //     sepset2[p] = sepset[k];
        //     p++;
        //   }
        //   bool sep_result = d_separated(level, n_node, i, j, sepset2, model);
        //   if (result && sep_result) {
        //     atomicAdd(stats, 1);
        //   } else if (result && !sep_result) {
        //     atomicAdd(stats + 1, 1);
        //   } else if (!result && sep_result) {
        //     atomicAdd(stats + 2, 1);
        //   } else {
        //     atomicAdd(stats + 3, 1);
        //   }
        // }
        __syncthreads();
      }
    }
  }
}

PDAG PCsearch(int citest_type, int n_node, int n_data,
              const vector<uint8_t> &data, const vector<int> &n_states,
              const vector<int> &model) {
  vector<int> G(n_node * n_node);
  for (int i = 0; i < n_node; i++) {
    for (int j = 0; j < n_node; j++) {
      if (i != j) G[i * n_node + j] = 1;
    }
  }
  vector<int> sepsets(n_node * n_node * max_level, -1);
  uint8_t *data_d;
  vector<double> regret(n_data * max_dim * 2);
  vector<int> stats(sm_num * 2);
  int *G_d, *n_states_d, *working_memory_d, *sepsets_d, *model_d, *stats_d;
  double *regret_d;
  int size_G = sizeof(int) * n_node * n_node;
  int size_data = sizeof(uint8_t) * n_data * n_node;
  int size_n_states = sizeof(int) * n_node;
  int size_working_memory = sizeof(int) * 500'000'000;
  int size_sepsets = sizeof(int) * n_node * n_node * max_level;
  int size_regret = sizeof(double) * n_data * max_dim * 2;
  int size_model = sizeof(int) * n_node * n_node * 2;
  int size_stats = sizeof(int) * sm_num * 2;
  CUDA_CHECK(cudaMalloc(&G_d, size_G));
  CUDA_CHECK(cudaMalloc(&data_d, size_data));
  CUDA_CHECK(cudaMalloc(&n_states_d, size_n_states));
  CUDA_CHECK(cudaMalloc(&working_memory_d, size_working_memory));
  CUDA_CHECK(cudaMalloc(&sepsets_d, size_sepsets));
  CUDA_CHECK(cudaMalloc(&regret_d, size_regret));
  CUDA_CHECK(cudaMalloc(&model_d, size_model));
  CUDA_CHECK(cudaMalloc(&stats_d, size_stats));
  CUDA_CHECK(
      cudaMemcpy(data_d, data.data(), size_data, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(n_states_d, n_states.data(), size_n_states,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(sepsets_d, sepsets.data(), size_sepsets,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(
      cudaMemcpy(model_d, model.data(), size_model, cudaMemcpyHostToDevice));
  CUDA_CHECK(
      cudaMemcpy(stats_d, stats.data(), size_stats, cudaMemcpyHostToDevice));

  calc_regret<<<n_data / 1024 + 1, 1024>>>(regret_d);
  CUDA_CHECK(
      cudaMemcpy(regret.data(), regret_d, size_regret, cudaMemcpyDeviceToHost));

  int true_max_deg = 0, true_max_indeg = 0, true_max_outdeg = 0;
  for (size_t i = 0; i < n_node; i++) {
    true_max_deg =
        max(true_max_deg,
            model.at(i * n_node) + model.at(n_node * n_node + i * n_node));
    true_max_indeg =
        max(true_max_indeg, model.at(n_node * n_node + i * n_node));
    true_max_outdeg = max(true_max_outdeg, model.at(i * n_node));
  }
  cout << "true_max_deg: " << true_max_deg
       << ", true_max_indeg: " << true_max_indeg
       << ", true_max_outdeg: " << true_max_outdeg << endl;
  // cout << "---regret---" << endl;
  // for (int n = n_data - 10; n <= n_data; n++) {
  //   for (int k = 1; k <= max_dim; k++) {
  //     cout << regret[n * max_dim + k - 1] << " ";
  //   }
  //   cout << endl;
  // }
  // cout << "---" << endl;
  vector<float> buf1;
  vector<int> buf2;
  cudaEvent_t start, stop;
  CUDA_CHECK(cudaEventCreate(&start));
  CUDA_CHECK(cudaEventCreate(&stop));
  // stage 1: Do CI tests between nodes and remove edge (undirected graph)
  int level = 0;
  int max_n_adj = n_node - 1;
  uint64_t max_dim_s = 1;
  while (level <= n_node - 2) {
    CUDA_CHECK(cudaMemcpy(G_d, G.data(), size_G, cudaMemcpyHostToDevice));
    stats = vector<int>(sm_num * 2);
    CUDA_CHECK(
        cudaMemcpy(stats_d, stats.data(), size_stats, cudaMemcpyHostToDevice));
    cout << "level: " << level << ", max_n_adj: " << max_n_adj << endl;
    if (level > max_level) break;
    CUDA_CHECK(cudaEventRecord(start));
    if (level == 0) {
      dim3 threadsPerBlock(64);
      dim3 numBlocks(n_node, n_node);
      PC_level_0<<<numBlocks, threadsPerBlock>>>(citest_type, n_node, n_data,
                                                 data_d, G_d, n_states_d,
                                                 regret_d, model_d, stats_d);
    } else {
      dim3 threadsPerBlock(64);
      dim3 numBlocks(n_node, max_n_adj * 2);
      uint64_t reserved_size_per_ci_test =
          max_dim_s * max_dim * max_dim + 2 * max_dim_s * max_dim + max_dim_s;
      if (reserved_size_per_ci_test > size_working_memory / sizeof(int)) {
        cout << "working memory is not enough" << endl;
        break;
      }
      if (reserved_size_per_ci_test * 2 < 1000) {
        PC_level_n<<<numBlocks, threadsPerBlock,
                     sizeof(int) * (max_n_adj + 2 + (level + 1) * 3 +
                                    reserved_size_per_ci_test) +
                         sizeof(double) * (5 + 1)>>>(
            citest_type, level, n_node, n_data, data_d, G_d, n_states_d, false,
            nullptr, sepsets_d, regret_d, model_d, stats_d);
      } else {
        uint64_t reserved_size_per_row =
            reserved_size_per_ci_test * numBlocks.y;
        int max_rows =
            size_working_memory / sizeof(int) / reserved_size_per_row;
        if (max_rows == 0) {
          int max_columns =
              size_working_memory / sizeof(int) / reserved_size_per_ci_test;
          numBlocks.x = 1;
          numBlocks.y = max_columns;
        } else if (numBlocks.x > max_rows) {
          numBlocks.x = max_rows;
        }
        cout << "numBlocks: " << numBlocks.x << ", " << numBlocks.y << endl;
        cout << "threadsPerBlock: " << threadsPerBlock.x << endl;
        PC_level_n<<<numBlocks, threadsPerBlock,
                     sizeof(int) * (max_n_adj + 2 + (level + 1) * 3) +
                         sizeof(double) * (5 + 1)>>>(
            citest_type, level, n_node, n_data, data_d, G_d, n_states_d, true,
            working_memory_d, sepsets_d, regret_d, model_d, stats_d);
      }
    }
    CUDA_CHECK(
        cudaMemcpy(stats.data(), stats_d, size_stats, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaEventRecord(stop));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));
    buf1.push_back(milliseconds / 1000);
    // cout << "----stats: " << endl;
    // for (int i = 0; i < sm_num; i++) {
    //   cout << i << ", " << stats[i] << ", " << stats[i + sm_num] << endl;
    // }
    // cout << "-----" << endl;
    buf2.push_back(accumulate(stats.begin(), stats.begin() + sm_num, 0));
    buf2.push_back(
        accumulate(stats.begin() + sm_num, stats.begin() + sm_num * 2, 0));
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
         << next_edge_cnt / 2 << endl;
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
  CUDA_CHECK(cudaEventDestroy(start));
  CUDA_CHECK(cudaEventDestroy(stop));
  for (float f : buf1) {
    cout << fixed << setprecision(3) << f << endl;
  }
  for (int c : buf2) {
    cout << c << endl;
  }
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

py::array_t<bool> gpuPC3(int citest_type, py::array_t<uint8_t> data,
                         py::array_t<int> n_states,
                         py::array_t<bool> true_model) {
  // translate input data to c++ vector(this is not optimal but I don't know
  // how to use pybind11::array_t)
  py::buffer_info buf_data = data.request(), buf_states = n_states.request(),
                  buf_model = true_model.request();
  const uint8_t *__restrict__ prt_data = static_cast<uint8_t *>(buf_data.ptr);
  const int *__restrict__ prt_states = static_cast<int *>(buf_states.ptr);
  const bool *__restrict__ prt_model = static_cast<bool *>(buf_model.ptr);
  size_t n_data = buf_data.shape[0],
         n_node = buf_data.shape[1];  // number of nodes
  cout << "n_data, n_node: " << n_data << ' ' << n_node << endl;
  vector<uint8_t> data_vec(n_data * n_node);
  vector<int> n_states_vec(n_node);
  vector<int> model_vec(2 * n_node * n_node);
  for (size_t i = 0; i < n_data; i++) {
    for (size_t j = 0; j < n_node; j++) {
      data_vec.at(j * n_data + i) = prt_data[i * n_node + j];
    }
  }
  for (size_t i = 0; i < n_node; i++) {
    n_states_vec.at(i) = prt_states[i];
  }
  for (size_t i = 0; i < n_node; i++) {
    int cnt = 0;
    for (size_t j = 0; j < n_node; j++) {
      if (prt_model[i * n_node + j]) {
        model_vec.at(i * n_node + (++cnt)) = j;
      }
    }
    model_vec.at(i * n_node) = cnt;
  }
  for (size_t i = 0; i < n_node; i++) {
    int cnt = 0;
    for (size_t j = 0; j < n_node; j++) {
      if (prt_model[j * n_node + i]) {
        model_vec.at(n_node * n_node + i * n_node + (++cnt)) = j;
      }
    }
    model_vec.at(n_node * n_node + i * n_node) = cnt;
  }
  auto endg = py::array_t<bool>({n_node, n_node});
  py::buffer_info buf_endg = endg.request();
  bool *__restrict__ prt_endg = static_cast<bool *>(buf_endg.ptr);

  PDAG Gend =
      PCsearch(citest_type, n_node, n_data, data_vec, n_states_vec, model_vec);

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