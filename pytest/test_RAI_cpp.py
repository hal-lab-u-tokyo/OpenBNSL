import pytest
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


@pytest.mark.parametrize(
    "data_type, sample_size, ess, max_iter, parallel, DP_threshold, search_neighbor, do_orientation_a2",
    [
        ("earthquake", 100000, -1, 1, 1, 1, True, True),  # 5 nodes
    ],
)

def test_RAI_cpp(
    data_type,
    sample_size,
    ess,
    max_iter,
    parallel,
    DP_threshold,
    search_neighbor,
    do_orientation_a2,
):
    calc_time = 0
    ave_score = [0, 0, 0, 0, 0, 0, 0]
    comparemodel = get_example_model(data_type)
    comparemodel = DAG2CPDAG(comparemodel)
    for i in range(max_iter):
        model = get_example_model(data_type)
        sampler = BayesianModelSampling(model)
        # sampler = GibbsSampling(model)
        data = sampler.forward_sample(size=sample_size, seed = 111)
        # data = model.simulate(sample_size)
        # print(data.head())
        best_model, raitime = RAIEstimator_cpp(
            data=data,
            ESS=ess,
            parallel=parallel,
            threshold_DP=DP_threshold,
            search_neighbor=search_neighbor,
            do_orientation_A2= do_orientation_a2,
        )
        best_model_compare = PDAG2CPDAG(best_model)
        calc_time += raitime
        # show(best_model_compare)
        errors = structural_errors(comparemodel, best_model_compare)
        print(
            f"iteration {i}: [SHD, ME, EE, DE, ED, MD, RD]:{errors[0]}, {errors[1]}, {errors[2]}, {errors[3]}, {errors[4]}, {errors[5]}, {errors[6]}"
        )
        ave_score = [x + y for x, y in zip(ave_score, errors)]
    ave_score = [x / max_iter for x in ave_score]
    calc_time /= max_iter
    print(
        f"[meanSHD, ME, EE, DE, ED, MD, RD]:{ave_score[0]}, {ave_score[1]}, {ave_score[2]}, {ave_score[3]}, {ave_score[4]}, {ave_score[5]}, {ave_score[6]}"
    )
    print(f"time: {calc_time}")

    if ave_score[0] > 0:
        print(f"Model: {data_type}, [SHD, ME, EE, DE, ED, MD, RD]: {ave_score}")
    assert ave_score[0] == 0, f"SHD too high for {data_type}: {ave_score[0]}"



