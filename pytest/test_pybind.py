import sys
import numpy as np
sys.path.append("/workspace/build")

import openbnsl
#dir(openbnsl)
c = openbnsl.myadd(1, 2)
print(c)

data=np.zeros((200,2))
for i in range(200):
    if i<50:
        data[i,0]=0
        data[i,1]=0
    elif i<100:
        data[i,0]=0
        data[i,1]=1
    elif i<150:
        data[i,0]=1
        data[i,1]=0
    else:
        data[i,0]=1
        data[i,1]=1

states=np.zeros(2)
states[0]=2
states[1]=2


endg = openbnsl.RAI(data, states, 1)
print(endg)