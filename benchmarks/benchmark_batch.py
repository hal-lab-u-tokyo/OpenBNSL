import subprocess

# models = ["diabetes", "munin"]
models = ["link"]
gpuPC_types = ["gpuPC", "gpuPC3"]
citest_types = ["g2", "sc"]
set_sizes = [20000, 40000, 100000, 200000, 400000]

for model in models:
  for gpuPC_type in gpuPC_types:
    for citest_type in citest_types:
      for set_size in set_sizes:
        subprocess.call("python3 -u benchmark_PC2.py "
        + gpuPC_type + " " + citest_type + " " + model + " " + str(set_size)
        + " |tee -a benchmark_result_link.txt", shell=True)