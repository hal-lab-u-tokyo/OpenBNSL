#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <iostream>
#include "template.h"
namespace py = pybind11;

py::array_t<double> matmul_naive(py::array_t<double> A, py::array_t<double> B) {
    py::buffer_info buf_A = A.request(), buf_B = B.request();
    if (buf_A.ndim != 2 || buf_B.ndim != 2)
        throw std::runtime_error("Number of dimensions must be two");
    if (buf_A.shape[1] != buf_B.shape[0])
        throw std::runtime_error("Input shapes must match");
    
    size_t M = buf_A.shape[0], K = buf_A.shape[1], N = buf_B.shape[1];
    auto C = py::array_t<double>({M, N});
    py::buffer_info buf_C = C.request();

    double *prt_A = static_cast<double *>(buf_A.ptr);
    double *prt_B = static_cast<double *>(buf_B.ptr);
    double *prt_C = static_cast<double *>(buf_C.ptr);

    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < N; j++) {
            double sum = 0;
            for (size_t k = 0; k < K; k++) {
                sum += prt_A[i * K + k] * prt_B[k * N + j];
            }
            prt_C[i * N + j] = sum;
        }
    }
    return C;
}

py::array_t<double> matmul_openmp(py::array_t<double> A, py::array_t<double> B) {
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

    #pragma omp parallel for
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < N; j++) {
            prt_C[i * N + j] = 0.0;
        }
    }

    #pragma omp parallel for 
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < K; j++) {
            const double a = prt_A[i * K + j];
            for (size_t k = 0; k < N; k++) {
                prt_C[i * N + k] += a * prt_B[j * N + k];
            }
        }
    }
    return C;
}
