import sys
sys.path.append("/workspace/build")

import openbnsl
#dir(openbnsl)
c = openbnsl.myadd(1, 2)
print(c)


G = openbnsl.RAI([1,1], 1.2)