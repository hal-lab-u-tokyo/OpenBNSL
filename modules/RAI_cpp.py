import sys
sys.path.append("/workspace/build")

import openbnsllib
import time
import numpy as np
from pgmpy.base import PDAG
sys.path.append('/workspace')
from modules.visualize_graph import display_graph_info as show





def RAIEstimator_cpp(data, ESS, parallel = 1, threshold_DP = 0, search_neighbor = True, do_orientation_A2 = True):

    #translate data from pandas dataframe to numpy array and extract column names
    columns = np.array(data.columns.tolist())
    t_data = data.to_numpy()
    #print(columns)

    #get statelist //重い
    statelist = [[] for i in range(len(columns))]
    for i in range(len(columns)):
        for j in range(t_data.shape[0]):
            if t_data[j][i] not in statelist[i]:
                statelist[i].append(t_data[j][i])

    #translate str matrix into int matrix and list of column names into list of int //重い
    data_int = np.zeros(t_data.shape, dtype = np.int8)
    for i in range(t_data.shape[0]):
        for j in range(len(columns)):
            data_int[i][j] = statelist[j].index(t_data[i][j])
    #print(data_int.dtype) 
    #print (t_data.dtype)
    #data_int = openbnsl.str2int_numpy(t_data, statelist)

    #テストコード
    #t_data = data.to_numpy(dtype=int)
    # for i in data.columns.tolist():
    #     k = 0
    #     for j in statelist[i]:
    #         data.replace({i: j}, int(k), inplace=True)
    #         k += 1
    # print(data)
    #data_int = data.to_numpy(dtype=int)
    #テストコード終了


    n_states = np.zeros(len(columns))
    for i in range(len(columns)):
        n_states[i] = len(statelist[i])
    #print(n_states)


    calc_time2 = 0
    t2 = time.time()
    ansmat = openbnsllib.structure_learning.RAI(data_int, n_states, ESS, parallel, threshold_DP, search_neighbor, do_orientation_A2)
    calc_time2 += time.time() - t2
    #print(f"shorttime: {calc_time2}")

    #print(ansmat)
    G = PDAG()
    G.add_nodes_from(columns)
    edges = []
    for i in range(ansmat.shape[0]):
        for j in range(ansmat.shape[1]):
            if ansmat[i][j] == 1:
                edges.append((columns[i], columns[j])) #i -> j
    G.add_edges_from(edges)
    #show(G)
    return G, calc_time2


