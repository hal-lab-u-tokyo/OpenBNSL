import timeit
import numpy as np
import pandas as pd
import os
import gc
import ctypes


class OpenMP:
    """
    Wrapper class for OpenMP functions using ctypes.
    Provides access to core OpenMP threading functionalities.
    """

    def __init__(self, library_name="libgomp.so.1"):
        try:
            self.lib = ctypes.CDLL(library_name)
        except OSError as e:
            raise OSError(f"Failed to load OpenMP library '{library_name}': {e}")

        # Define function prototypes
        self.lib.omp_get_max_threads.restype = ctypes.c_int
        self.lib.omp_get_num_threads.restype = ctypes.c_int
        self.lib.omp_set_num_threads.argtypes = [ctypes.c_int]

    def get_max_threads(self):
        """Returns the maximum number of threads available."""
        return self.lib.omp_get_max_threads()

    def get_num_threads(self):
        """Returns the number of threads in the current parallel region."""
        return self.lib.omp_get_num_threads()

    def set_num_threads(self, num_threads):
        """Sets the number of threads to use in parallel regions."""
        if not isinstance(num_threads, int) or num_threads <= 0:
            raise ValueError("Number of threads must be a positive integer.")
        self.lib.omp_set_num_threads(num_threads)


def benchmark_functions_variants(
    func_variants, args=None, trials=10, unit="ms", benchmark_name="Benchmark"
):
    """
    Benchmarks multiple variants of the same function with the same input over a number of trials.

    Parameters:
        func_variants (list of tuples): List of (variant_name, function_callable, thread_count) triples.
        args (tuple): Arguments to pass to all function variants.
        trials (int): Number of times to run each function for benchmarking.
        unit (str): Time unit for output (options: "s", "ms", "us", "ns").
        benchmark_name (str): Name of the benchmark to display in the output.

    Returns:
        None: Prints the benchmarking results in a table format.
    """
    # Time conversion factors
    unit_factors = {
        "s": 1,
        "ms": 1000,
        "us": 1_000_000,
        "ns": 1_000_000_000,
    }

    if unit not in unit_factors:
        raise ValueError(
            f"Invalid unit '{unit}'. Supported units: {', '.join(unit_factors.keys())}."
        )

    factor = unit_factors[unit]
    results = []

    # Benchmark each function variant
    omp = OpenMP()
    max_threads = omp.get_max_threads()
    for variant_name, func, thread_count in func_variants:
        if thread_count > max_threads:
            raise ValueError(
                f"Thread count {thread_count} > max threads {max_threads}."
            )
        if thread_count > 0:
            omp.set_num_threads(thread_count)
        else:
            omp.set_num_threads(1)

        times = []
        for _ in range(trials):
            gc.collect()  # Clear cache before each trial
            time_taken = timeit.timeit(lambda: func(*args), number=1)
            times.append(time_taken)

        times = np.array(times)
        times_converted = times * factor  # Convert times to the desired unit

        # Calculate statistics
        q1, q3 = np.percentile(times, [25, 75])
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        lower_outliers = np.sum(times < lower_bound)
        upper_outliers = np.sum(times > upper_bound)

        min_time = np.min(times_converted)
        max_time = np.max(times_converted)
        mean_time = np.mean(times_converted)
        stddev_time = np.std(times_converted)
        median_time = np.median(times_converted)
        iqr_time = (q3 - q1) * factor
        ops = factor / mean_time if mean_time > 0 else float("inf")

        # Append results for this function variant
        results.append(
            {
                "Variant (time in {})".format(unit): variant_name,
                "Threads": thread_count if thread_count > 0 else "N/A",
                "Min": f"{min_time:.2f}",
                "Max": f"{max_time:.2f}",
                "Mean": f"{mean_time:.2f}",
                "StdDev": f"{stddev_time:.2f}",
                "Median": f"{median_time:.2f}",
                "IQR": f"{iqr_time:.2f}",
                "Outliers": f"{lower_outliers};{upper_outliers}",
                "OPS": f"{ops:.2f}",
                "Rounds": trials,
                "Iterations": 1,
            }
        )

    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)

    # Calculate dynamic header width
    table_width = len(results_df.to_string(index=False, header=True).split("\n")[0])
    benchmark_header = f" {benchmark_name}: {len(func_variants)} variants "
    header_width = max(len(benchmark_header), table_width)
    header_line = "-" * header_width

    # Display the table
    print(header_line)
    print(benchmark_header.center(header_width))
    print(header_line)
    print(results_df.to_string(index=False, header=True))
    print(header_line)
