import sys

import openbnsllib
import time
import numpy as np
from pgmpy.base import PDAG

from modules.visualize_graph import display_graph_info as show
from modules.dataframe_wrapper2 import dataframe_wrapper2


def RAIEstimator_cpp(
    data, ESS, parallel=1, threshold_DP=0, search_neighbor=True, do_orientation_A2=True
):
    columns, data_int, n_states = dataframe_wrapper2(data)

    calc_time2 = 0
    t2 = time.perf_counter()
    ansmat = openbnsllib.structure_learning.RAI(
        data_int,
        n_states,
        ESS,
        parallel,
        threshold_DP,
        search_neighbor,
        do_orientation_A2,
    )
    calc_time2 += time.perf_counter() - t2
    # print(f"shorttime: {calc_time2}")
    # print(ansmat)
    G = PDAG()
    G.add_nodes_from(columns)
    edges = []
    for i in range(ansmat.shape[0]):
        for j in range(ansmat.shape[1]):
            if ansmat[i][j] == 1:
                edges.append((columns[i], columns[j]))  # i -> j
    G.add_edges_from(edges)
    # show(G)
    return G, calc_time2
