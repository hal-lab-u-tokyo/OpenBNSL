import pytest
import random
import csv
import os
import time
import pandas as pd
from pgmpy.utils import get_example_model
import logging

logger = logging.getLogger(__name__)

import openbnsllib
from helpers.omp import OpenMP
from helpers.data_cache import get_samples_head
from helpers.pgmpy_bridge import to_pgmpy, to_openbnsl
from helpers.structural_distance import structural_errors

RESULTS_PATH = os.path.join("benchmarks", "results", "pc_and_rai")
SCENARIO_NAME = os.path.splitext(os.path.basename(__file__))[0]
RESULTS_FILE = os.path.join(RESULTS_PATH, f"{SCENARIO_NAME}.csv")
SUMMARY_FILE = os.path.join(RESULTS_PATH, f"{SCENARIO_NAME}_summary.csv")


def initialize():
    if not os.path.exists(RESULTS_PATH):
        os.makedirs(RESULTS_PATH)
    with open(RESULTS_FILE, "w", newline="") as f:
        csv.writer(f).writerow(
            [
                "model_name",
                "algo",
                "num_threads",
                "num_samples",
                "seed",
                "shd",
                "original_score",
                "learned_score",
                "score_ratio",
                "elapsed_sec",
            ]
        )


def summarize():
    df = pd.read_csv(RESULTS_FILE)
    summary = (
        df.groupby(["model_name", "algo", "num_threads", "num_samples"])
        .agg(
            count=("shd", "count"),
            shd_mean=("shd", "mean"),
            shd_std=("shd", "std"),
            original_score_mean=("original_score", "mean"),
            # original_score_std=("original_score", "std"),
            learned_score_mean=("learned_score", "mean"),
            # learned_score_std=("learned_score", "std"),
            score_ratio_mean=("score_ratio", "mean"),
            # score_ratio_std=("score_ratio", "std"),
            time_mean_s=("elapsed_sec", "mean"),
            time_std_s=("elapsed_sec", "std"),
            # time_min_s=("elapsed_sec", "min"),
            # time_max_s=("elapsed_sec", "max"),
        )
        .reset_index()
        # .sort_values(["model_name", "algo", "num_threads"])
    )
    with pd.option_context("display.float_format", "{:.2f}".format):
        print(f"\n--- Summary for {SCENARIO_NAME} ---")
        print(summary.to_string(index=False))
    summary.to_csv(SUMMARY_FILE, index=False)


ALGORITHMS = {
    "pc": lambda dfw, ci, cols: openbnsllib.structure_learning.pc(
        dfw, ci, max_cond_vars=len(cols)
    ),
    "rai": lambda dfw, ci, cols: openbnsllib.structure_learning.rai(
        dfw, ci, max_cond_vars=len(cols)
    ),
}


@pytest.mark.parametrize(
    "model_name",
    [
        "alarm",  # test
        # "asia", "cancer", "earthquake", "sachs", "survey",  # Small networks
        # "alarm", "barley", "child", "insurance", "mildew", "water",  # Medium networks
        # "hailfinder", "hepar2", "win95pts", # Large networks
        # "andes", "diabetes", "link", "munin1", "pathfinder", "pigs", # Extra large networks
        # "munin", "munin2", "munin3", "munin4", # Very extra large networks
    ],
)
@pytest.mark.parametrize("algo", list(ALGORITHMS.keys()))
@pytest.mark.parametrize("num_threads", [1, 16])
# @pytest.mark.parametrize("num_samples", [int(1e4), int(2e5), int(2e6)])
@pytest.mark.parametrize("num_samples", [int(1e4)])
@pytest.mark.parametrize("seed", [0])
def benchmark_compare_algos(model_name, algo, num_threads, num_samples, seed):

    # Setup
    random.seed(seed)
    omp = OpenMP()
    omp.set_num_threads(num_threads)

    original_pdag_pgmpy = get_example_model(model_name)

    def _gen(num_samples_: int, seed_: int):
        return original_pdag_pgmpy.simulate(num_samples_, seed=seed_)

    samples = get_samples_head(
        model_name=model_name,
        seed=seed,
        num_samples=num_samples,
        generator=_gen,
    )

    cols = list(samples.columns)
    df_wrapper = openbnsllib.base.DataframeWrapper(samples)
    citest_type = openbnsllib.citest.ChiSquare(level=0.01)
    original_pdag_obnsl = to_openbnsl(original_pdag_pgmpy, df_wrapper.col_str2idx)

    # Trial
    start = time.perf_counter()
    learned_pdag_obnsl = ALGORITHMS[algo](df_wrapper, citest_type, cols)
    elapsed = time.perf_counter() - start

    # Measurement
    learned_pdag_pgmpy = to_pgmpy(learned_pdag_obnsl, cols)
    error_dict = structural_errors(original_pdag_pgmpy, learned_pdag_pgmpy)
    shd = error_dict["SHD"]

    original_score = original_pdag_obnsl.score(df_wrapper, openbnsllib.score.BDeu(1.0))
    learned_score = learned_pdag_obnsl.score(df_wrapper, openbnsllib.score.BDeu(1.0))
    score_ratio = (
        learned_score / original_score if original_score != 0 else float("inf")
    )

    print(
        f"[{algo.upper()}] model={model_name}, num_threads={num_threads}, "
        f"num_samples={num_samples}, seed={seed}, "
        f"shd={shd}, original_score={original_score:.2f}, "
        f"learned_score={learned_score:.2f}, "
        f"score_ratio={score_ratio:.2f}, elapsed={elapsed:.2f}s"
    )

    with open(RESULTS_FILE, "a", newline="") as f:
        csv.writer(f).writerow(
            [
                model_name,
                algo,
                num_threads,
                num_samples,
                seed,
                shd,
                original_score,
                learned_score,
                score_ratio,
                elapsed,
            ]
        )

    assert True
