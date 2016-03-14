import numpy as np
import numpy.linalg as LA
from inputGF import *
from surfGF import *

g = Input()


#klist = [[0,0.041667,0]]

eng_list = [4.1414, 4.4914, 4.6419]


surf_spectral=[]
surf_modified_spectral = []

for kpt in g.klist:
    H00, H01 = construct_H00_H01(g, kpt)

    epsilon0 = H00
    alpha0 = H01
    beta0 = H01.conj().T
    tmp, tmp1 = per_k(g,kpt)
    surf_spectral += [tmp]
    surf_modified_spectral += [tmp1]
surf_spectral = np.array(surf_spectral)
surf_modified_spectral = np.array(surf_modified_spectral)
print(surf_spectral, surf_modified_spectral)

