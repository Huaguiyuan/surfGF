import numpy as np
import numpy.linalg as LA
from input import *
from surfGF import *

g = Input()
#print(g.b)
kpt = [0,0.041667,0]
#kpt = [0.1 ,0.2 ,0]
#print(g.kpoint_scale(kpt,g.b))
#print(kpoint_scale(kpt,g.b))
#print(g.eng_list)

#print(g.dict, g.a,g.b)

H00, H01 = construct_H00_H01(g, kpt)

epsilon0 = H00
alpha0 = H01
beta0 = H01.conj().T

#print(H00,'\n',H01)
'''
print(np.imag(H00))
x=0
for i in range(60):
   x+=np.imag(H00)[i,i]
print(x)
'''
#print(LA.eigvals(H00))

#print(max(abs(H00.flatten())))
#print(max(abs(H01.flatten())))
#eng_list = [1.40510429]
eng_list = [4.1414, 4.4914, 4.6419]
surf_spectralN=[]

for i in range(len(eng_list)):

    omega = eng_list[i]
    epsilon, epsilons = iteration.iteration.iterate(alpha0, beta0, epsilon0, epsilon0,omega,smearing = g.dict['smearing'], prec = g.dict['convergence'], max_step=g.dict['maximum_iteration'])
    surf_spectralN.append(spectral_weight(g,omega,epsilons, smearing = g.dict['smearing']))

print(surf_spectralN)

