# benchmarks/scenarios/benchmark_template.py
import pytest
import random
import csv
import os
import time
import pandas as pd
from pgmpy.utils import get_example_model
import openbnsllib
from helpers.omp import OpenMP
from helpers.pgmpy_bridge import to_pgmpy
from helpers.structural_distance import structural_errors

RESULTS_PATH = os.path.join("benchmarks", "results")
SCENARIO_NAME = os.path.splitext(os.path.basename(__file__))[0]
RESULTS_FILE = os.path.join(RESULTS_PATH, f"{SCENARIO_NAME}.csv")
SUMMARY_FILE = os.path.join(RESULTS_PATH, f"{SCENARIO_NAME}_summary.csv")


def initialize():
    """
    Prepare a fresh results CSV for this scenario.
    """
    if not os.path.exists(RESULTS_PATH):
        os.makedirs(RESULTS_PATH)
    with open(RESULTS_FILE, "w", newline="") as f:
        csv.writer(f).writerow(
            ["model_name", "num_threads", "seed", "shd", "elapsed_sec"]
        )


@pytest.mark.parametrize("model_name", ["cancer", "asia", "child", "alarm"])
@pytest.mark.parametrize("num_threads", [1, 2, 4, 8])
@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def benchmark_template(model_name, num_threads, seed):
    """
    One benchmark trial.
    """

    # Setup
    random.seed(seed)
    omp = OpenMP()
    omp.set_num_threads(num_threads)

    original_pdag = get_example_model(model_name)
    samples = original_pdag.simulate(int(1e3), seed=seed)
    samples = samples[sorted(samples.columns)]
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    citest_type = openbnsllib.citest.ChiSquare(0.01)

    # Trial
    start = time.perf_counter()
    learned_pdag = openbnsllib.structure_learning.pc(
        df_wrapper, citest_type, max_cond_vars=len(samples.columns), stable=True
    )
    elapsed = time.perf_counter() - start

    # Measurement
    error_dict = structural_errors(
        original_pdag, to_pgmpy(learned_pdag, list(samples.columns))
    )
    shd = error_dict["SHD"]

    with open(RESULTS_FILE, "a", newline="") as f:
        csv.writer(f).writerow([model_name, num_threads, seed, shd, elapsed])

    assert True


def summarize():
    df = pd.read_csv(RESULTS_FILE)
    summary = (
        df.groupby(["model_name", "num_threads"])
        .agg(
            count=("shd", "count"),
            shd_mean=("shd", "mean"),
            shd_std=("shd", "std"),
            time_mean_s=("elapsed_sec", "mean"),
            time_std_s=("elapsed_sec", "std"),
            time_min_s=("elapsed_sec", "min"),
            time_max_s=("elapsed_sec", "max"),
        )
        .reset_index()
    )
    print(f"\n--- Summary for {SCENARIO_NAME} ---")
    print(summary.to_string(index=False))

    summary.to_csv(SUMMARY_FILE, index=False)
