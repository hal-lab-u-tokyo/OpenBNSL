import pandas as pd
import numpy as np
import networkx as nx
import argparse
import os
from pgmpy.utils import get_example_model

from pgmpy.sampling import BayesianModelSampling
from pgmpy.estimators import BicScore, HillClimbSearch
from pgmpy.estimators import PC

import time
import sys 
sys.path.append('/workspace')
from modules.RAIEstimator_fixed import RAIEstimator, RAIEstimator_transitivity
from modules.structural_distance import structural_errors, DAG2CPDAG, PDAG2CPDAG
from modules.visualize_graph import display_graph_info as show
from modules.CITests_fixed import NatoriScore
from build.CMakeFiles.CMakeTm
SAVE_DIR = "./results_2"

ESTIMATOR={
    "RAI": RAIEstimator,
    "RAI_t": RAIEstimator_transitivity,
    "HC": HillClimbSearch,
    "PC": PC
}
SCORE={
    "BIC": BicScore
}

def load_data(type, sample_size):
    model = get_example_model(type)
    sampler = BayesianModelSampling(model)
    data = sampler.forward_sample(size=sample_size)
    #print(data.head())
    return model, data

def test_benchmark(
        estimate_type,
        data_type,
        sample_size,
        structure_score,
        ess,
        max_iter,
    ): 
    calc_time = 0
    ave_score = [0, 0, 0, 0, 0, 0, 0] 
    comparemodel = get_example_model(data_type)
    comparemodel = DAG2CPDAG(comparemodel)
    for i in range(max_iter):
        model , data = load_data(data_type, sample_size)
        print(data.head())
        data = 


def save_benchmark(estimate_type, data_type, size, structure_score, ave_score, calc_time):
    results = pd.DataFrame(columns=["estimate_type", "data_type", "structure_score", "ave_score", "calc_time"])
    results.loc[0] = [estimate_type, data_type, structure_score, ave_score, calc_time]
    if not os.path.exists(SAVE_DIR):
        os.makedirs(SAVE_DIR)
    file_path = f"{SAVE_DIR}/{estimate_type}_{data_type}_{size}_{structure_score}"
    file_exists = os.path.isfile(file_path)
    with open(file_path, "w") as f:
        results.to_csv(f, header=not file_exists, index=False)


def arg_parser():
    parser = argparse.ArgumentParser(description="Benchmarking for Bayesian Network Structure Learning")
    parser.add_argument("--estimate_type", type=str, default="RAI", help="Estimator type")
    parser.add_argument("--data_type", type=str, default="andes", help="Data type")
    parser.add_argument("--sample_size", type=int, default=10000, help="Sample size")
    parser.add_argument("--structure_score", type=str, default="BIC", help="Structure score")
    parser.add_argument("--ess", type=float, default=5, help="ESS")
    parser.add_argument("--max_iter", type=int, default=1, help="Number of iterations")
    return parser.parse_args()

def main():
    args = arg_parser()
    data_type = args.data_type
    sample_size = args.sample_size
    max_iter = args.max_iter
    estimate_type = args.estimate_type
    ess = args.ess
    structure_score = args.structure_score
    file_path = f"{SAVE_DIR}/{estimate_type}_{data_type}_{sample_size}_{structure_score}"
    test_benchmark(estimate_type,
                   data_type,
                   sample_size,
                   structure_score,
                   ess,
                   max_iter)



if __name__ == "__main__":
    main()
