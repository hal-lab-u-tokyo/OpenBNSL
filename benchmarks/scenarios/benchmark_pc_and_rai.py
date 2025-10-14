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

RESULTS_PATH = os.path.join(
    "benchmarks", "results", "pc_and_rai", time.strftime("%Y%m%d-%H%M%S")
)
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
                "num_vars",
                "num_samples",
                "citest",
                "algo",
                "num_threads",
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
        df.groupby(
            ["model_name", "num_vars", "num_samples", "citest", "algo", "num_threads"]
        )
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
        .sort_values(
            ["model_name", "num_vars", "num_samples", "citest", "algo", "num_threads"]
        )
    )
    with pd.option_context("display.float_format", "{:.2f}".format):
        print(f"\n--- Summary for {SCENARIO_NAME} ---")
        print(summary.to_string(index=False))
    summary.to_csv(SUMMARY_FILE, index=False)


ALGORITHMS = {
    "pc_edge_parallel": lambda dfw, ci, cols, timeout_sec: openbnsllib.structure_learning.pc_edge_parallel(
        dfw, ci, max_cond_vars=len(cols), timeout_sec=timeout_sec
    ),
    "rai_edge_parallel": lambda dfw, ci, cols, timeout_sec: openbnsllib.structure_learning.rai_edge_parallel(
        dfw, ci, max_cond_vars=len(cols), timeout_sec=timeout_sec
    ),
}
CITESTS = {
    "chi2": openbnsllib.citest.ChiSquare,
    # "g2": openbnsllib.citest.GSquare,
}


@pytest.mark.parametrize(
    "model_name",
    [
        "asia",
        "cancer",
        "earthquake",
        "sachs",
        "survey",
        "alarm",
        # "child", "insurance", "water", "hailfinder", "win95pts",
        # "barley", "mildew", "hepar2", "andes", "munin1", "diabetes",
        # "link", "munin", "munin2", "munin3", "munin4",
        # "pathfinder", "pigs",
    ],
)
@pytest.mark.parametrize("citest", list(CITESTS.keys()))
@pytest.mark.parametrize("algo", list(ALGORITHMS.keys()))
# @pytest.mark.parametrize("num_threads", [1, 16, 128])
@pytest.mark.parametrize("num_threads", [128])
# @pytest.mark.parametrize("num_samples", [int(1e4), int(1e5)])
@pytest.mark.parametrize("num_samples", [int(1e4)])
@pytest.mark.parametrize("timeout_sec", [3600])
# @pytest.mark.parametrize("seed", [0,1,2,3,4])
@pytest.mark.parametrize("seed", [0])
def benchmark_compare_algos(
    model_name, citest, algo, num_threads, num_samples, timeout_sec, seed
):

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
    num_vars = df_wrapper.num_vars
    # oracle_graph = to_openbnsl(original_pdag_pgmpy, df_wrapper.col_str2idx)
    # citest_type = openbnsllib.citest.OracleGraph(oracle_graph)
    citest_type = CITESTS[citest](0.05)
    original_pdag_obnsl = to_openbnsl(original_pdag_pgmpy, df_wrapper.col_str2idx)

    # Trial
    print(
        f"model={model_name}, num_vars={num_vars}, num_samples={num_samples}, algo={algo.upper()}, citest={citest.upper()}, num_threads={num_threads}, seed={seed} ..."
    )
    start = time.perf_counter()
    learned_pdag_obnsl = ALGORITHMS[algo](
        df_wrapper, citest_type, cols, timeout_sec=timeout_sec
    )
    elapsed = time.perf_counter() - start

    # Measurement
    print("Evaluating ...")
    learned_pdag_pgmpy = to_pgmpy(learned_pdag_obnsl, cols)
    error_dict = structural_errors(original_pdag_pgmpy, learned_pdag_pgmpy)
    shd = error_dict["SHD"]

    original_score = original_pdag_obnsl.score(df_wrapper, openbnsllib.score.BDeu(1.0))
    learned_score = learned_pdag_obnsl.score(df_wrapper, openbnsllib.score.BDeu(1.0))
    score_ratio = (
        learned_score / original_score if original_score != 0 else float("inf")
    )

    print(error_dict)

    print(
        f"model={model_name}, num_vars={num_vars}, num_samples={num_samples}, citest={citest.upper()}, algo={algo.upper()}, num_threads={num_threads}, seed={seed}, "
        f"shd={shd}/{num_vars * (num_vars - 1) // 2} (worst case), "
        f"original_score={original_score:.2f}, "
        f"learned_score={learned_score:.2f}, "
        f"score_ratio={score_ratio:.2f}, elapsed={elapsed:.2f}s"
    )

    with open(RESULTS_FILE, "a", newline="") as f:
        csv.writer(f).writerow(
            [
                model_name,
                num_vars,
                num_samples,
                citest,
                algo,
                num_threads,
                seed,
                shd,
                original_score,
                learned_score,
                score_ratio,
                elapsed,
            ]
        )

    assert True
