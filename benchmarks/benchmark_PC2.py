from pgmpy.utils import get_example_model
from pgmpy.base import PDAG

import sys, time
import numpy as np

sys.path.append("/workspace")
sys.path.append("/workspace/build")

import openbnsllib
from modules.structural_distance import structural_errors, PDAG2CPDAG, retrieve_adjacency_matrix


def get_list(line):
  return list(map(int, line[:-1].split(",")))

def main():
    path = "/workspace/dataset/"
    model_name = "munin"
    set_size = 20000
    set_id = 5
    file_name = path + model_name + "_" + str(set_id) + "_int_s400000"
    n_states = []
    data_int = []

    with open(file_name, mode='r') as f:
      lines = f.readlines()
      n_states = np.array(get_list(lines[0]))
      for i in range(1, set_size + 1):
        data_int.append(get_list(lines[i]))
      data_int = np.array(data_int)
    
    model = get_example_model(model_name)
    model_mat = retrieve_adjacency_matrix(model).astype(dtype=bool)
    nodes = list(model.nodes)

    print("timer start")
    start_time = time.perf_counter()
    ansmat = openbnsllib.structure_learning.gpuPC3(data_int, n_states, model_mat)
    calc_time = time.perf_counter() - start_time
    print("ansmat = ", ansmat)
    print("calc_time = ", calc_time)

    G = PDAG()
    G.add_nodes_from(nodes)
    edges = []
    for i in range(ansmat.shape[0]):
        for j in range(ansmat.shape[1]):
            if ansmat[i][j] == 1:
                edges.append((nodes[i], nodes[j]))  # i -> j
    G.add_edges_from(edges)
    errors = structural_errors(
        PDAG2CPDAG(model), PDAG2CPDAG(G)
    )
    errors.append(calc_time)
    print("errors:", errors)
    print(model, G)


if __name__ == "__main__":
    main()
