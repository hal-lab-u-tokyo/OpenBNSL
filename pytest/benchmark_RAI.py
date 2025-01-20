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

sys.path.append("/workspace")
from modules.RAIEstimator_fixed import RAIEstimator, RAIEstimator_transitivity
from modules.RAI_cpp import RAIEstimator_cpp
from modules.PC_cpp import PCEstimator_cpp
from modules.structural_distance import structural_errors, DAG2CPDAG, PDAG2CPDAG
from modules.visualize_graph import display_graph_info as show
from modules.CITests_fixed import NatoriScore

SAVE_DIR = "./save_benchmark_RAI.txt"


def load_data(type, sample_size):
    model = get_example_model(type)
    sampler = BayesianModelSampling(model)
    # sampler = GibbsSampling(model)
    data = sampler.forward_sample(size=sample_size)
    # data = model.simulate(sample_size)
    # print(data.head())
    return model, data


def do_benchmark(
    data_type,
    sample_size,
    max_iter,
    parallel,
    DP,
    search_neighbor,
    do_orientation_A2,
    f,
):
    calc_time = 0
    ave_score = [0, 0, 0, 0, 0, 0, 0]
    comparemodel = get_example_model(data_type)
    comparemodel = DAG2CPDAG(comparemodel)
    DP_threshold = 0
    if DP == 1:
        DP_threshold = 15
    for i in range(max_iter):
        model, data = load_data(data_type, sample_size)
        best_model, shorttime = RAIEstimator_cpp(
            data=data,
            ESS=-10,
            parallel=parallel,
            threshold_DP=DP_threshold,
            search_neighbor=search_neighbor,
            do_orientation_A2=do_orientation_A2,
        )
        calc_time += shorttime
        best_model_compare = PDAG2CPDAG(best_model)
        # show(best_model_compare)
        errors = structural_errors(comparemodel, best_model_compare)
        ave_score = [x + y for x, y in zip(ave_score, errors)]
    ave_score = [x / max_iter for x in ave_score]
    calc_time /= max_iter
    # print(
    #     f"[meanSHD, ME, EE, DE, ED, MD, RD]:{ave_score[0]}, {ave_score[1]}, {ave_score[2]}, {ave_score[3]}, {ave_score[4]}, {ave_score[5]}, {ave_score[6]}"
    # )
    # print(f"time: {calc_time}")
    f.write(
        f"[{data_type},{sample_size},{parallel},{DP},{search_neighbor},{do_orientation_A2}][meanSHD, ME, EE, DE, ED, MD, RD]:{ave_score[0]}, {ave_score[1]}, {ave_score[2]}, {ave_score[3]}, {ave_score[4]}, {ave_score[5]}, {ave_score[6]}"
    )
    f.write(f", time: {calc_time}\n")
    return

def compare_algorithm(data_type, sample_size, max_iter, f):
    parallel = 0
    DP = 0
    search_neighbor = False
    do_orientation_A2 = False
    do_benchmark(data_type, sample_size, max_iter, parallel, DP, search_neighbor, do_orientation_A2, f) # C++ RAI
    do_orientation_A2 = True
    do_benchmark(data_type, sample_size, max_iter, parallel, DP, search_neighbor, do_orientation_A2, f) # C++ RAI + orientation A2
    do_orientation_A2 = False
    search_neighbor = True
    do_benchmark(data_type, sample_size, max_iter, parallel, DP, search_neighbor, do_orientation_A2, f) # C++ RAI + search neighbor
    do_orientation_A2 = True
    do_benchmark(data_type, sample_size, max_iter, parallel, DP, search_neighbor, do_orientation_A2, f) # C++ RAI + orientation A2 + search neighbor
    DP = 1
    do_benchmark(data_type, sample_size, max_iter, parallel, DP, search_neighbor, do_orientation_A2, f) # C++ RAI + orientation A2 + search neighbor + DP
    DP = 0
    parallel = 1
    do_benchmark(data_type, sample_size, max_iter, parallel, DP, search_neighbor, do_orientation_A2, f) # C++ RAI + orientation A2 + search neighbor + parallel
    DP = 1
    do_benchmark(data_type, sample_size, max_iter, parallel, DP, search_neighbor, do_orientation_A2, f) # C++ RAI + orientation A2 + search neighbor + DP + parallel
    return



def macroinstruction(f):
    max_iter = 1
    data_types = [
        "cancer",
        "earthquake",
        "survey",
        "sachs",
        "child",
        "alarm",
        #  "win95pts",
    ]  # "win95"はやらない
    sample_sizes = [10000, 20000, 50000, 100000, 200000, 1000000, 2000000, 10000000]
    parallel = 0
    DP = 0
    search_neighbor = False
    do_orientation_A2 = True
    f.write("\n")
    for data_type in data_types:
        for sample_size in sample_sizes:
            if (
                (data_type == "cancer" and sample_size == 50000)
                or (data_type == "cancer" and sample_size == 100000)
                or (data_type == "alarm" and sample_size == 50000)
                or (data_type == "alarm" and sample_size == 100000)
                or (data_type == "alarm" and sample_size == 20000)
                or (data_type == "alarm" and sample_size == 1000000)
                or (data_type == "alarm" and sample_size == 10000000)
                or (data_type == "survey" and sample_size == 2000000)
                or (data_type == "sachs" and sample_size == 20000)
                or (data_type == "sachs" and sample_size == 100000)
                or (data_type == "sachs" and sample_size == 2000000)
                or (data_type == "child" and sample_size == 1000000)
                or (data_type == "child" and sample_size == 10000000)
                or (data_type == "win95pts" and sample_size == 20000)
                or (data_type == "win95pts" and sample_size == 100000)
                or (data_type == "win95pts" and sample_size == 1000000)
                or (data_type == "win95pts" and sample_size == 2000000)
                or (data_type == "win95pts" and sample_size == 10000000)
            ):
                continue
            compare_algorithm(data_type, sample_size, max_iter, f)


    return


def main():
    f = open(SAVE_DIR, "w")
    f.write("start benchmark\n")
    macroinstruction(f)
    f.close()


if __name__ == "__main__":
    main()
