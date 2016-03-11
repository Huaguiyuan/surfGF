import numpy as np
import numpy.linalg as LA
from read_input import *
from surfGF import *

g = GF()
surf_spectral = []

g.klist = [[0,0,0],[1,1,1]]
for kpt in g.klist : #g.dict['klist']:
    surf_spectral += [per_k(g,kpt)]

##supercell problem

###write result

print(len(g.klist), len(g.eng_list), len(surf_spectral), len(surf_spectral[0]))
file = open('surface.dat', 'w')
for i in range(len(g.klist)):
    for j in range(len(g.eng_list)):
        print('{0:4f}    {1:8f}    {2:8f}    {3:8f}    {4:8f}    {5:8f}    {6:8f}'.format(i, g.kd[i], g.klist[i][0], g.klist[i][1], g.klist[i][2], g.eng_list[j], surf_spectral[i][j]), file=file)
file.close()


