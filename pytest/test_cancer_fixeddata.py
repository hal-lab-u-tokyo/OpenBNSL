import sys
import numpy as np
sys.path.append("/workspace/build")

import openbnsl
#dir(openbnsl)

#test data from cancer n=1000000, [Pollution, Smoke, Cancer, X-ray]
data=np.zeros((1000000,4))
count = 0
n_data = [503496, 125874, 63, 567, 209520, 52380, 810, 7290, 54880, 13720, 140, 1260, 22800, 5700, 150, 1350]
for i in range(n_data[0]):
    data[count,0]=0
    data[count,1]=0
    data[count,2]=0
    data[count,3]=0
    count += 1
for i in range(n_data[1]):
    data[count,0]=0
    data[count,1]=0
    data[count,2]=0
    data[count,3]=1
    count += 1
for i in range(n_data[2]):
    data[count,0]=0
    data[count,1]=0
    data[count,2]=1
    data[count,3]=0
    count += 1
for i in range(n_data[3]):
    data[count,0]=0
    data[count,1]=0
    data[count,2]=1
    data[count,3]=1
    count += 1
for i in range(n_data[4]):
    data[count,0]=0
    data[count,1]=1
    data[count,2]=0
    data[count,3]=0
    count += 1
for i in range(n_data[5]):
    data[count,0]=0
    data[count,1]=1
    data[count,2]=0
    data[count,3]=1
    count += 1
for i in range(n_data[6]):
    data[count,0]=0
    data[count,1]=1
    data[count,2]=1
    data[count,3]=0
    count += 1
for i in range(n_data[7]):
    data[count,0]=0
    data[count,1]=1
    data[count,2]=1
    data[count,3]=1
    count += 1
for i in range(n_data[8]):
    data[count,0]=1
    data[count,1]=0
    data[count,2]=0
    data[count,3]=0
    count += 1
for i in range(n_data[9]):
    data[count,0]=1
    data[count,1]=0
    data[count,2]=0
    data[count,3]=1
    count += 1
for i in range(n_data[10]):
    data[count,0]=1
    data[count,1]=0
    data[count,2]=1
    data[count,3]=0
    count += 1
for i in range(n_data[11]):
    data[count,0]=1
    data[count,1]=0
    data[count,2]=1
    data[count,3]=1
    count += 1
for i in range(n_data[12]):
    data[count,0]=1
    data[count,1]=1
    data[count,2]=0
    data[count,3]=0
    count += 1
for i in range(n_data[13]):
    data[count,0]=1
    data[count,1]=1
    data[count,2]=0
    data[count,3]=1
    count += 1
for i in range(n_data[14]):
    data[count,0]=1
    data[count,1]=1
    data[count,2]=1
    data[count,3]=0
    count += 1
for i in range(n_data[15]):
    data[count,0]=1
    data[count,1]=1
    data[count,2]=1
    data[count,3]=1
    count += 1


states=np.ones(4)
states[0]=2
states[1]=2
states[3]=2
states[2]=2

endg = openbnsl.RAI(data, states)
print(endg)

