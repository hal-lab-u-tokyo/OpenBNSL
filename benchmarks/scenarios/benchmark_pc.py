import pytest
import random
import csv
import os
import time
import pandas as pd
from pgmpy.utils import get_example_model
from pgmpy.estimators import PC as PgmpyPC
import logging

logger = logging.getLogger(__name__)

import openbnsllib
from helpers.omp import OpenMP
from helpers.data_cache import get_samples_head
from helpers.pgmpy_bridge import to_pgmpy, to_openbnsl
from helpers.structural_distance import structural_errors

RESULTS_PATH = os.path.join(
    "benchmarks", "results", "pc_vs_pgmpy", time.strftime("%Y%m%d-%H%M%S")
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
            learned_score_mean=("learned_score", "mean"),
            score_ratio_mean=("score_ratio", "mean"),
            time_mean_s=("elapsed_sec", "mean"),
            time_std_s=("elapsed_sec", "std"),
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


ALGORITHMS = ["pc_edge_parallel", "pc_pgmpy"]

CITESTS = {
    "chi2": openbnsllib.citest.ChiSquare,
}

PGMPY_CITEST_MAP = {
    "chi2": "chi_square",
}


@pytest.mark.parametrize(
    "model_name",
    [
        "alarm",
        # "barley",
        "child",
        "insurance",
        # "mildew",
        "water",
    ],
)
@pytest.mark.parametrize("citest", list(CITESTS.keys()))
@pytest.mark.parametrize("algo", ALGORITHMS)
@pytest.mark.parametrize("num_threads", [128])
@pytest.mark.parametrize("num_samples", [int(2e5)])
@pytest.mark.parametrize("timeout_sec", [3600])
@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
def benchmark_compare_algos(
    model_name, citest, algo, num_threads, num_samples, timeout_sec, seed
):

    random.seed(seed)
    omp = OpenMP()
    omp.set_num_threads(num_threads)

    original_pdag_pgmpy = get_example_model(model_name)

    def _gen(num_samples_, seed_):
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

    citest_type = CITESTS[citest](0.05)

    original_pdag_obnsl = to_openbnsl(original_pdag_pgmpy, df_wrapper.col_str2idx)

    print(
        f"model={model_name}, num_vars={num_vars}, num_samples={num_samples}, "
        f"algo={algo.upper()}, citest={citest.upper()}, num_threads={num_threads}, seed={seed} ..."
    )

    start = time.perf_counter()

    if algo == "pc_edge_parallel":
        learned_pdag_obnsl = openbnsllib.structure_learning.pc_edge_parallel(
            df_wrapper,
            citest_type,
            max_cond_vars=len(cols),
            timeout_sec=timeout_sec,
        )
        learned_pdag_pgmpy = to_pgmpy(learned_pdag_obnsl, cols)

    elif algo == "pc_pgmpy":
        pgmpy_ci_name = PGMPY_CITEST_MAP[citest]

        est = PgmpyPC(samples)
        learned_pdag_pgmpy = est.estimate(
            variant="parallel",
            ci_test=pgmpy_ci_name,
            return_type="pdag",
            significance_level=0.05,
            max_cond_vars=len(cols),
            n_jobs=num_threads,
            show_progress=False,
        )

        learned_pdag_obnsl = to_openbnsl(learned_pdag_pgmpy, df_wrapper.col_str2idx)

    else:
        raise ValueError(f"Unknown algo: {algo}")

    elapsed = time.perf_counter() - start

    print("Evaluating ...")

    error_dict = structural_errors(original_pdag_pgmpy, learned_pdag_pgmpy)
    shd = error_dict["SHD"]

    original_score = original_pdag_obnsl.score(df_wrapper, openbnsllib.score.BDeu(1.0))
    learned_score = learned_pdag_obnsl.score(df_wrapper, openbnsllib.score.BDeu(1.0))

    score_ratio = (
        learned_score / original_score if original_score != 0 else float("inf")
    )

    print(error_dict)
    print(
        f"model={model_name}, num_vars={num_vars}, num_samples={num_samples}, "
        f"citest={citest.upper()}, algo={algo.upper()}, num_threads={num_threads}, seed={seed}, "
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
