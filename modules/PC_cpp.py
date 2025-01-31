import sys
sys.path.append("/workspace/build")

import openbnsl
import numpy as np
from pgmpy.base import PDAG
sys.path.append('/workspace')
from modules.visualize_graph import display_graph_info as show





def PCEstimator_cpp(data, ESS):
    #translate data from pandas dataframe to numpy array and extract column names
    columns = np.array(data.columns.tolist())
    t_data = data.to_numpy()
    #print(columns)

    #get statelist
    statelist = [[] for i in range(len(columns))]
    for i in range(len(columns)):
        for j in range(t_data.shape[0]):
            if t_data[j][i] not in statelist[i]:
                statelist[i].append(t_data[j][i])
    #print(statelist)

    #translate str matrix into int matrix and list of column names into list of int
    data_int = np.zeros(t_data.shape)
    for i in range(t_data.shape[0]):
        for j in range(len(columns)):
            data_int[i][j] = statelist[j].index(t_data[i][j])
    #print(data_int)

    n_states = np.zeros(len(columns))
    for i in range(len(columns)):
        n_states[i] = len(statelist[i])
    #print(n_states)

    ansmat = openbnsl.PC(data_int, n_states, ESS)

    G = PDAG()
    G.add_nodes_from(columns)
    edges = []
    for i in range(ansmat.shape[0]):
        for j in range(ansmat.shape[1]):
            if ansmat[i][j] == 1:
                edges.append((columns[i], columns[j])) #i -> j
    G.add_edges_from(edges)
    #show(G)
    return G


