import numpy as np
DTe_Dxq = np.zeros([12,13])
for i in range(0,6):
    DTe_Dxq[i, i] = 1
for i in range(9,12):
    DTe_Dxq[i,i+1] = 1
print(DTe_Dxq)