from pgmpy.readwrite import BIFReader
from pgmpy.base import PDAG

import sys, time, pickle
import numpy as np

sys.path.append("/workspace")
sys.path.append("/workspace/build")

import openbnsllib
from modules.structural_distance import structural_errors, PDAG2CPDAG, retrieve_adjacency_matrix


def get_list(line):
  return list(map(int, line[:-1].split(",")))

def main():
    path = "/workspace/dataset/"
    gpuPC_type = sys.argv[1]
    citest_type = 0 if sys.argv[2] == "g2" else 1
    # model_name = "diabetes"
    # set_size = 20000
    model_name = sys.argv[3]
    set_size = int(sys.argv[4])
    print(gpuPC_type, citest_type, model_name, set_size)
    for set_id in range(10):
      print("set_id: " + str(set_id))
      file_name = path + model_name + "_" + str(set_id) + "_int_s400000"
      data_int = []

      with open(file_name, mode='r') as f:
        lines = f.readlines()
        n_states = np.array(get_list(lines[0]))
        for i in range(1, set_size + 1):
          data_int.append(get_list(lines[i]))
        data_int = np.array(data_int)
      
      # reader = BIFReader(path + model_name + ".bif")
      # model = reader.get_model()
      with open(path + model_name + ".pickle", "rb") as f:
        model = pickle.load(f)
      model_mat = retrieve_adjacency_matrix(model).astype(dtype=bool)
      nodes = list(model.nodes)

      print("timer start")
      start_time = time.perf_counter()
      if gpuPC_type == "gpuPC":
        ansmat = openbnsllib.structure_learning.gpuPC(citest_type, data_int, n_states, model_mat)
      else:
        ansmat = openbnsllib.structure_learning.gpuPC3(citest_type, data_int, n_states, model_mat)
      calc_time = time.perf_counter() - start_time
      # print("ansmat = ", ansmat)
      # print("calc_time = ", calc_time)

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
      errors.append(round(calc_time, 3))
      for e in errors:
        print(e)
      # print(model, G)


if __name__ == "__main__":
    main()
