#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <vector>

#include "template.h"
namespace py = pybind11;

__global__ void matmul(const double* A,
                       const double* B,
                       double* C,
                       int M,
                       int K,
                       int N) {
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  int col = blockIdx.x * blockDim.x + threadIdx.x;

  if (row < M && col < N) {
    double val = 0.0;
    for (int i = 0; i < K; ++i) {
      val += A[row * K + i] * B[i * N + col];
    }
    C[row * N + col] = val;
  }
}

#define CUDA_CHECK(call)                               \
  do {                                                 \
    cudaError_t e = call;                              \
    if (e != cudaSuccess) {                            \
      throw std::runtime_error(cudaGetErrorString(e)); \
    }                                                  \
  } while (0)

py::array_t<double> matmul_cuda(py::array_t<double> A, py::array_t<double> B) {
  py::buffer_info buf_A = A.request(), buf_B = B.request();
  if (buf_A.ndim != 2 || buf_B.ndim != 2)
    throw std::runtime_error("Number of dimensions must be two");
  if (buf_A.shape[1] != buf_B.shape[0])
    throw std::runtime_error("Input shapes must match");

  size_t M = buf_A.shape[0], K = buf_A.shape[1], N = buf_B.shape[1];
  auto C = py::array_t<double>({M, N});
  py::buffer_info buf_C = C.request();

  const double* __restrict__ prt_A = static_cast<double*>(buf_A.ptr);
  const double* __restrict__ prt_B = static_cast<double*>(buf_B.ptr);
  double* __restrict__ prt_C = static_cast<double*>(buf_C.ptr);

  // Allocate device memory
  double *d_A, *d_B, *d_C;
  CUDA_CHECK(cudaMalloc(&d_A, M * K * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_B, K * N * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_C, M * N * sizeof(double)));

  // Copy data to device
  CUDA_CHECK(
      cudaMemcpy(d_A, prt_A, M * K * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(
      cudaMemcpy(d_B, prt_B, K * N * sizeof(double), cudaMemcpyHostToDevice));

  // Launch kernel
  dim3 threads(16, 16);
  dim3 blocks((N + threads.x - 1) / threads.x, (M + threads.y - 1) / threads.y);
  matmul<<<blocks, threads>>>(d_A, d_B, d_C, M, K, N);

  // wait for kernel to finish
  CUDA_CHECK(cudaDeviceSynchronize());

  // Copy result back to host
  CUDA_CHECK(
      cudaMemcpy(prt_C, d_C, M * N * sizeof(double), cudaMemcpyDeviceToHost));

  // Free device memory
  CUDA_CHECK(cudaFree(d_A));
  CUDA_CHECK(cudaFree(d_B));
  CUDA_CHECK(cudaFree(d_C));

  return C;
}