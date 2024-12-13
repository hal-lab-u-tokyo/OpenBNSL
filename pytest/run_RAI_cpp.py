import sys
sys.path.append("/workspace/build")

import openbnsl
# c = openbnsl.addadd(1, 2)
# print(c)


import pandas as pd
import numpy as np
import networkx as nx
import argparse
import os
from pgmpy.utils import get_example_model
from pgmpy.base import PDAG
from pgmpy.sampling import BayesianModelSampling
from pgmpy.estimators import BicScore, HillClimbSearch
from pgmpy.estimators import PC

import time
import sys 
sys.path.append('/workspace')
from modules.RAIEstimator_fixed import RAIEstimator, RAIEstimator_transitivity
#from modules.RAI_cpp import RAIEstimator_cpp
from modules.structural_distance import structural_errors, DAG2CPDAG, PDAG2CPDAG
from modules.visualize_graph import display_graph_info as show
from modules.CITests_fixed import NatoriScore
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

def translate_data(data):
    columns = data.columns.tolist()
    data = data.values.tolist()
    #print(data)
    return data, columns

def RAIEstimator_cpp(data, ESS):
    t_data, columns = translate_data(data)
    ansmat = openbnsl.RAI(t_data, ESS)
    G = PDAG()
    G.add_nodes_from(columns)
    edges = []
    for i in range(ansmat.shape[0]):
        for j in range(ansmat.shape[1]):
            if ansmat[i][j] == 1:
                edges.append((columns[i], columns[j]))
    G.add_edges_from(edges)
    return G

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
        data2, columns = translate_data(data)
        estimator = ESTIMATOR[estimate_type](data)
        score = SCORE[structure_score](data)
        t = time.time()
        best_model = RAIEstimator_cpp(data, ESS=ess)
        show(best_model)
    #     #best_model = estimator.estimate(ESS = ess)
    #     calc_time += time.time() - t
    #     #compare_best_model = PDAG2CPDAG(best_model)
    #     errors = structural_errors(comparemodel, best_model)
    #     print(f"iteration {i}: [SHD, ME, EE, DE, ED, MD, RD]:{errors[0]}, {errors[1]}, {errors[2]}, {errors[3]}, {errors[4]}, {errors[5]}, {errors[6]}")
    #     #print(f"iteration {i}: SHD:{errors[0]}, ME:{errors[1]}, EE:{errors[2]}, {errors[3]}, {errors[4]}, {errors[5]}, {errors[6]}")
    #     ave_score = [x + y for x, y in zip(ave_score, errors)]
    # ave_score = [x / max_iter for x in ave_score]
    # calc_time /= max_iter
    # print(f"[meanSHD, ME, EE, DE, ED, MD, RD]:{ave_score[0]}, {ave_score[1]}, {ave_score[2]}, {ave_score[3]}, {ave_score[4]}, {ave_score[5]}, {ave_score[6]}")
    # print(f"time: {calc_time}")
    # save_benchmark(estimate_type, data_type, sample_size, structure_score, ave_score, calc_time)


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
    parser.add_argument("--estimate_type", type=str, default="RAI_t", help="Estimator type")
    parser.add_argument("--data_type", type=str, default="cancer", help="Data type")
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
