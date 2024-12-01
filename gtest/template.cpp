#include "template.h"

#include <gtest/gtest.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <random>
#include <vector>

namespace py = pybind11;

TEST(MatmulTest, correctness) {
  const int n = 3, m = 4, p = 5;
  const double epsilon = 1e-6;
  const int seed = 42;

  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(0, 1);

  auto a = py::array_t<double>({n, m});
  auto b = py::array_t<double>({m, p});
  auto a_buf = a.request();
  auto b_buf = b.request();
  double *a_ptr = static_cast<double *>(a_buf.ptr);
  double *b_ptr = static_cast<double *>(b_buf.ptr);

  for (int i = 0; i < n * m; i++) a_ptr[i] = dist(gen);
  for (int i = 0; i < m * p; i++) b_ptr[i] = dist(gen);

  auto c_naive = matmul_naive(a, b);
  auto c_openmp = matmul_openmp(a, b);
  auto c_naive_buf = c_naive.request();
  auto c_openmp_buf = c_openmp.request();
  double *c_naive_ptr = static_cast<double *>(c_naive_buf.ptr);
  double *c_openmp_ptr = static_cast<double *>(c_openmp_buf.ptr);

  for (int i = 0; i < n * p; i++) {
    ASSERT_NEAR(c_naive_ptr[i], c_openmp_ptr[i], epsilon);
  }
}
