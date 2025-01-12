import sys
import numpy as np
import math
sys.path.append("/workspace/build")

import openbnsl
#dir(openbnsl)

#test data from cancer n=1000000, [Pollution, Smoke, Cancer, X-ray]
#data=np.zeros((1000000,4))
data=np.zeros((2000000,4))
count = 0
n_data = [440559*2, 188811*2, 441, 819, 183330*2, 78570*2, 2835*2, 5265*2, 48020*2, 20580*2, 490*2, 910*2, 19950*2, 8550*2, 525*2, 975*2]
print(n_data)
#n_data = [503496, 125874, 63, 567, 209520, 52380, 810, 7290, 54880, 13720, 140, 1260, 22800, 5700, 150, 1350]
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
ESS = 10
endg = openbnsl.RAI(data, states, ESS, 1, 10)
print(endg)

dependent_score = 0
independent_score = 0
alpha = 0
# dependent_score =  dependent_score + 2 * math.lgamma(alpha * 4) - math.lgamma(alpha * 4 + 1391859)- math.lgamma(alpha * 4 + 608141) - 8 * math.lgamma(alpha)
# dependent_score =  dependent_score + math.lgamma(alpha + 881559)+ math.lgamma(alpha + 378441)+ math.lgamma(alpha + 372330)+ math.lgamma(alpha + 167670)+ math.lgamma(alpha + 97020)+ math.lgamma(alpha + 42980)+ math.lgamma(alpha + 40950)+ math.lgamma(alpha + 19050)

# independent_score =  independent_score + 4 * math.lgamma(alpha * 2) - math.lgamma(alpha * 2 + 1391859)- math.lgamma(alpha * 2 + 608141) - math.lgamma(alpha * 2 + 1391859)- math.lgamma(alpha * 2 + 608141) - 8 * math.lgamma(alpha)
# independent_score =  independent_score + math.lgamma(alpha + 1253889) + math.lgamma(alpha + 546111) + math.lgamma(alpha + 137970) + math.lgamma(alpha + 978579) + math.lgamma(alpha + 421421) + math.lgamma(alpha + 62030) + math.lgamma(alpha + 186720) + math.lgamma(alpha + 413280) 

# dependent_score =  dependent_score - math.lgamma(1391859)- math.lgamma(608141)
# dependent_score =  dependent_score + math.lgamma(alpha + 881559)+ math.lgamma(alpha + 378441)+ math.lgamma(372330)+ math.lgamma(alpha + 167670)+ math.lgamma(alpha + 97020)+ math.lgamma(alpha + 42980)+ math.lgamma(alpha + 40950)+ math.lgamma(alpha + 19050)

# independent_score =  independent_score - math.lgamma(alpha * 2 + 1391859)- math.lgamma(alpha * 2 + 608141) - math.lgamma(alpha * 2 + 1391859)- math.lgamma(alpha * 2 + 608141)
# independent_score =  independent_score + math.lgamma(alpha + 1253889) + math.lgamma(alpha + 546111) + math.lgamma(alpha + 137970) + math.lgamma(alpha + 978579) + math.lgamma(alpha + 421421) + math.lgamma(alpha + 62030) + math.lgamma(alpha + 186720) + math.lgamma(alpha + 413280) 


print(f"dependent_score: {dependent_score}")
print(f"independent_score: {independent_score}")