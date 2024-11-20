import os
import sys
import time
import numpy as np

sys.path.append("/workspace/build")
import openbn

def measure_time(func, *args, trials=5):
    times = []
    for _ in range(trials):
        start = time.time()
        func(*args)
        end = time.time()
        times.append(end - start)
    return np.mean(times), np.std(times)

def python_matmul(a, b):
    n, m = len(a), len(b[0])
    p = len(b)
    result = [[0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for k in range(p):
                result[i][j] += a[i][k] * b[k][j]
    return result

trials = 10
matrix_size = 256
np.random.seed(42)
a = np.random.rand(matrix_size, matrix_size)
b = np.random.rand(matrix_size, matrix_size)

print(f"Matrix size: {matrix_size}x{matrix_size}")
print(f"Number of trials: {trials}")

python_time, python_std = measure_time(python_matmul, a.tolist(), b.tolist(), trials=trials)
print(f"Python Time: {python_time:.6f} ± {python_std:.6f} seconds")

naive_time, naive_std = measure_time(openbn.matmul_naive, a, b, trials=trials)
print(f"OpenBN Naive Time: {naive_time:.6f} ± {naive_std:.6f} seconds")

openmp_time, openmp_std = measure_time(openbn.matmul_openmp, a, b, trials=trials)
print(f"OpenBN OpenMP Time: {openmp_time:.6f} ± {openmp_std:.6f} seconds")

numpy_time, numpy_std = measure_time(np.dot, a, b, trials=trials)
print(f"NumPy Time: {numpy_time:.6f} ± {numpy_std:.6f} seconds")


