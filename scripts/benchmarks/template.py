import sys
import numpy as np
from benchmark_utils import benchmark_functions_variants

sys.path.append("/workspace/build")
import openbnsl


def python_matmul(a, b):
    n, m = len(a), len(b[0])
    p = len(b)
    result = [[0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for k in range(p):
                result[i][j] += a[i][k] * b[k][j]
    return result


matrix_size = 256
np.random.seed(0)
a = np.random.rand(matrix_size, matrix_size)
b = np.random.rand(matrix_size, matrix_size)
input_args = (a, b)
func_variants = [
    ("Python", python_matmul, -1),
    ("Naive C++", openbnsl.matmul_naive, -1),
    ("OpenMP", openbnsl.matmul_openmp, 1),
    ("OpenMP", openbnsl.matmul_openmp, 4),
    ("OpenMP", openbnsl.matmul_openmp, 8),
    ("NumPy", np.dot, -1),
]
if openbnsl.cuda_enabled:
    func_variants += [("CUDA", openbnsl.matmul_cuda, -1)]

benchmark_functions_variants(
    func_variants,
    input_args,
    trials=30,
    unit="ms",
    benchmark_name="Benchmark Matrix Multiplication ({}x{})".format(a.shape, b.shape),
)
