from pgmpy.utils import get_example_model
from pgmpy.sampling import BayesianModelSampling
import numpy as np

model_name = "diabetes"
sample_size = 400000
set_cnt = 10

model = get_example_model(model_name)
for set_id in range(set_cnt):
  path = "/workspace/dataset/"

  sampler = BayesianModelSampling(model)
  data = sampler.forward_sample(sample_size)
  t_data = data.to_numpy()
  columns = np.array(data.columns.tolist())

  file_name = path + model_name + "_" + str(set_id) +  "_s" + str(sample_size)
  with open(file_name, mode='w') as f:
    line = ",".join(list(columns))
    f.write(line + "\n")
    for i in range(t_data.shape[0]):
      line = ",".join(list(t_data[i, :]))
      f.write(line + "\n")

  statelist = [[] for i in range(len(columns))]

  for i in range(len(columns)):
    for j in range(t_data.shape[0]):
      if t_data[j][i] not in statelist[i]:
        statelist[i].append(t_data[j][i])
  for i in range(len(columns)):
    statelist[i].sort()
  print("max_dim: {}".format(max([len(x) for x in statelist])))

  data_int = np.zeros(t_data.shape, dtype=np.uint8)
  for i in range(t_data.shape[0]):
    for j in range(len(columns)):
      data_int[i][j] = statelist[j].index(t_data[i][j])

  n_states = np.zeros(len(columns), dtype=np.int32)
  for i in range(len(columns)):
    n_states[i] = len(statelist[i])

  file_name2 = path + model_name + "_" + str(set_id) + "_int_s" + str(sample_size)
  with open(file_name2, mode='w') as f:
    line = ",".join(map(str, list(n_states)))
    f.write(line + "\n")
    for i in range(t_data.shape[0]):
      line = ",".join(map(str, list(data_int[i, :])))
      f.write(line + "\n")