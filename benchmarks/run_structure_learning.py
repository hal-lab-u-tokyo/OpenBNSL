import pandas as pd
import numpy as np
import networkx as nx
import argparse
import os
from pgmpy.utils import get_example_model

from pgmpy.sampling import BayesianModelSampling, GibbsSampling
from pgmpy.estimators import BicScore, HillClimbSearch
from pgmpy.estimators import PC

import time
import sys

from modules.RAIEstimator_fixed import RAIEstimator, RAIEstimator_transitivity
from modules.RAI_cpp import RAIEstimator_cpp
from modules.PC_cpp import PCEstimator_cpp
from modules.structural_distance import _structural_errors, DAG2CPDAG, PDAG2CPDAG
from modules.CITests_fixed import NatoriScore

SAVE_DIR = "./results_2"

ESTIMATOR = {
    "RAI": RAIEstimator,
    "RAI_cpp": RAIEstimator_cpp,
    "PC_cpp": PCEstimator_cpp,
    "RAI_t": RAIEstimator_transitivity,
    "HC": HillClimbSearch,
    "PC": PC,
}
SCORE = {"BIC": BicScore}


def load_data(type, sample_size):
    model = get_example_model(type)
    sampler = BayesianModelSampling(model)
    # sampler = GibbsSampling(model)
    data = sampler.forward_sample(size=sample_size)
    # data = model.simulate(sample_size)
    # print(data.head())
    return model, data


def test_benchmark(
    estimate_type,
    data_type,
    sample_size,
    structure_score,
    ess,
    max_iter,
    parallel,
    DP_threshold,
    search_neighbor,
):
    calc_time = 0
    ave_score = [0, 0, 0, 0, 0, 0, 0]
    comparemodel = get_example_model(data_type)
    comparemodel = DAG2CPDAG(comparemodel)
    for i in range(max_iter):
        model, data = load_data(data_type, sample_size)
        if not (estimate_type == "RAI_cpp" or estimate_type == "PC_cpp"):
            estimator = ESTIMATOR[estimate_type](data)
        score = SCORE[structure_score](data)
        t = time.time()
        if estimate_type == "RAI_cpp":
            best_model, shorttime = RAIEstimator_cpp(
                data=data,
                ESS=ess,
                parallel=parallel,
                threshold_DP=DP_threshold,
                search_neighbor=search_neighbor,
            )
        elif estimate_type == "PC_cpp":
            best_model = PCEstimator_cpp(data=data, ESS=ess)
        elif estimate_type == "RAI":
            best_model = estimator.estimate(ESS=ess)
        elif estimate_type == "RAI_t":
            best_model = estimator.estimate(ESS=ess)
        elif estimate_type == "HC":
            best_model = estimator.estimate(scoring_method=score)
        elif estimate_type == "PC":
            best_model = estimator.estimate()
        else:
            raise ValueError("Invalid estimator type")
        if estimate_type == "RAI_cpp":
            calc_time += shorttime
        else:
            calc_time += time.time() - t
        best_model_compare = PDAG2CPDAG(best_model)
        errors = _structural_errors(comparemodel, best_model_compare)
        print(
            f"iteration {i}: [SHD, ME, EE, DE, ED, MD, RD]:{errors[0]}, {errors[1]}, {errors[2]}, {errors[3]}, {errors[4]}, {errors[5]}, {errors[6]}"
        )
        # print(f"iteration {i}: SHD:{errors[0]}, ME:{errors[1]}, EE:{errors[2]}, {errors[3]}, {errors[4]}, {errors[5]}, {errors[6]}")
        ave_score = [x + y for x, y in zip(ave_score, errors)]
    ave_score = [x / max_iter for x in ave_score]
    calc_time /= max_iter
    print(
        f"[meanSHD, ME, EE, DE, ED, MD, RD]:{ave_score[0]}, {ave_score[1]}, {ave_score[2]}, {ave_score[3]}, {ave_score[4]}, {ave_score[5]}, {ave_score[6]}"
    )
    print(f"time: {calc_time}")
    save_benchmark(
        estimate_type, data_type, sample_size, structure_score, ave_score, calc_time
    )


def save_benchmark(
    estimate_type, data_type, size, structure_score, ave_score, calc_time
):
    results = pd.DataFrame(
        columns=[
            "estimate_type",
            "data_type",
            "structure_score",
            "ave_score",
            "calc_time",
        ]
    )
    results.loc[0] = [estimate_type, data_type, structure_score, ave_score, calc_time]
    if not os.path.exists(SAVE_DIR):
        os.makedirs(SAVE_DIR)
    file_path = f"{SAVE_DIR}/{estimate_type}_{data_type}_{size}_{structure_score}"
    file_exists = os.path.isfile(file_path)
    with open(file_path, "w") as f:
        results.to_csv(f, header=not file_exists, index=False)


def arg_parser():
    parser = argparse.ArgumentParser(
        description="Benchmarking for Bayesian Network Structure Learning"
    )
    parser.add_argument(
        "--estimate_type", type=str, default="RAI_cpp", help="Estimator type"
    )
    parser.add_argument("--data_type", type=str, default="alarm", help="Data type")
    parser.add_argument("--sample_size", type=int, default=20000, help="Sample size")
    parser.add_argument(
        "--structure_score", type=str, default="BIC", help="Structure score"
    )
    parser.add_argument("--ess", type=float, default=5, help="ESS")
    parser.add_argument("--max_iter", type=int, default=1, help="Number of iterations")
    parser.add_argument("--parallel", type=int, default=1, help="parallel")
    parser.add_argument("--DP_threshold", type=int, default=15, help="DP_threshold")
    parser.add_argument(
        "--search_neighbor", type=bool, default=False, help="search_neighbor"
    )
    return parser.parse_args()


def main():
    args = arg_parser()
    data_type = args.data_type
    sample_size = args.sample_size
    max_iter = args.max_iter
    estimate_type = args.estimate_type
    ess = args.ess
    parallel = args.parallel
    DP_threshold = args.DP_threshold
    structure_score = args.structure_score
    search_neighbor = args.search_neighbor
    file_path = (
        f"{SAVE_DIR}/{estimate_type}_{data_type}_{sample_size}_{structure_score}"
    )
    test_benchmark(
        estimate_type,
        data_type,
        sample_size,
        structure_score,
        ess,
        max_iter,
        parallel,
        DP_threshold,
        search_neighbor,
    )


if __name__ == "__main__":
    main()
