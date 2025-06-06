while True:
  test_cond = input().split()
  if len(test_cond) == 0:
    break
  if test_cond[0] == "gpuPC":
    metrics_len = 16
  else:
    metrics_len = 20
  lst = [[] for i in range(metrics_len)]
  for i in range(10):
    assert int(input().split()[1]) == i
    line = input()
    while ":" in line or line == "timer start":
      line = input()
    lst[0].append(float(line))
    for j in range(1, metrics_len):
      lst[j].append(float(input()))
  print(test_cond)
  for i in range(metrics_len):
    print(",".join(map(str, lst[i])))
