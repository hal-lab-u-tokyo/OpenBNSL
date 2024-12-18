import sys
import numpy as np
sys.path.append("/workspace/build")

import openbnsl
#dir(openbnsl)


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


endg = openbnsl.RAI(data, states)
print(endg)


# import sys
# import numpy as np
# sys.path.append("/workspace/build")

# import openbnsl
# #dir(openbnsl)


# data=np.zeros((40,3))
# for i in range(40):
#     if i<5:
#         data[i,0]=0
#         data[i,1]=0
#         data[i,2]=0
#     elif i<10:
#         data[i,0]=0
#         data[i,1]=0
#         data[i,2]=0
#     elif i<15:
#         data[i,0]=0
#         data[i,1]=1
#         data[i,2]=0
#     elif i<20:
#         data[i,0]=0
#         data[i,1]=1
#         data[i,2]=1
#     elif i<25:
#         data[i,0]=1
#         data[i,1]=0
#         data[i,2]=0
#     elif i<30:
#         data[i,0]=1
#         data[i,1]=0
#         data[i,2]=1
#     elif i<35:
#         data[i,0]=1
#         data[i,1]=1
#         data[i,2]=1
#     else:
#         data[i,0]=1
#         data[i,1]=1
#         data[i,2]=1

# states=np.zeros(3)
# states[0]=2
# states[1]=2
# states[2]=2

# endg = openbnsl.RAI(data, states, 1)
# print(endg)