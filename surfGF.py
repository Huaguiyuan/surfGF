import numpy as np
import numpy.linalg as LA
from read_input import *


def point_scale(pt, a):
    # point_scale, including rpt and kpt, if it is rpt, put a as lattice, if it is kpt, put a as inversed lattice
    pt_scaled = a[0, :] * pt[0] + a[1, :] * pt[1] + a[2, :] * pt[2]
    pt_scaled = np.array(pt_scaled, dtype='float')
    return pt_scaled


def construct_H00_H01(g, kpt):
    principle_layer_thickness = g.dict['principle_layer_thickness']
    surface_index = g.dict['slab_direction']

    plane_dist = []  # plane distance to the plane which crosses (0, 0, 0)
    for i in range(g.nrpt):
        plane_dist.append(
            surface_index[0] * g.rpt[i, 0] + surface_index[1] * g.rpt[i, 1] + surface_index[2] * g.rpt[i, 2])

    max_dist = max(plane_dist)
    min_dist = min(plane_dist)

    atom_layer = []
    atom_layer1 = []
    for i in range(2 * principle_layer_thickness - 1):
        j = i - principle_layer_thickness + 1
        if j <= max_dist and j >= min_dist:
            atom_layer.append([k for k, x in enumerate(plane_dist) if x == j])
        else:
            atom_layer.append('None')
        if i + 1 <= max_dist:
            atom_layer1.append([k for k, x in enumerate(plane_dist) if x == i + 1])
        else:
            atom_layer1.append('None')
    '''
    print(g.rpt)
    print(atom_layer)
    print(atom_layer1)
    '''

    H00 = np.zeros((principle_layer_thickness * g.num_wann, principle_layer_thickness * g.num_wann), dtype=np.complex)
    H01 = np.zeros((principle_layer_thickness * g.num_wann, principle_layer_thickness * g.num_wann), dtype=np.complex)
    for i in range(principle_layer_thickness):
        for j in range(principle_layer_thickness):
            hamij = np.zeros((g.num_wann, g.num_wann), dtype=np.complex)
            hamij1 = np.zeros((g.num_wann, g.num_wann), dtype=np.complex)
            num_atom_layer = j - i + principle_layer_thickness - 1  # number between j and i
            num_atom_layer1 = num_atom_layer

            if atom_layer[num_atom_layer] != 'None':
                for k in atom_layer[num_atom_layer]:
                    hamij += np.exp(1j * np.dot(point_scale(kpt, g.b), point_scale(g.rpt[k], g.a))) * (
                    g.hamr[:, :, k] + 1j * g.hami[:, :, k])  # *g.weight[k]
            H00[g.num_wann * i:g.num_wann * (i + 1), g.num_wann * j:g.num_wann * (j + 1)] = hamij

            if atom_layer1[num_atom_layer1] != 'None':
                for k1 in atom_layer1[num_atom_layer1]:
                    hamij1 += np.exp(1j * np.dot(point_scale(kpt, g.b), point_scale(g.rpt[k1], g.a))) * (
                    g.hamr[:, :, k1] + 1j * g.hami[:, :, k1])  # *g.weight[k]
            H01[g.num_wann * i:g.num_wann * (i + 1), g.num_wann * j:g.num_wann * (j + 1)] = hamij1

    return H00, H01


def iterative_method(alpha0, beta0, epsilon0, epsilons0, omega, smearing=0.01, precision=0.01, min_step=12, max_step=20,
                     flag=0):
    m, n = np.shape(epsilon0)
    Omega = np.eye(m, dtype=np.complex) * (omega + smearing * 1j)
    green = LA.inv(Omega - epsilon0)
    alpha = np.dot(alpha0, np.dot(green, alpha0))
    beta = np.dot(beta0, np.dot(green, beta0))
    epsilon = epsilon0 + np.dot(alpha0, np.dot(green, beta0)) + np.dot(beta0, np.dot(green, alpha0))
    epsilons = epsilons0 + np.dot(alpha0, np.dot(green, beta0))

    while np.max(np.abs(epsilon - epsilon0)) > precision and np.max(
            np.abs(epsilons - epsilons0)) > precision and flag < max_step or flag < min_step:
        flag += 1
        return iterative_method(alpha, beta, epsilon, epsilons, omega, smearing, precision, min_step, max_step, flag)
    return epsilon, epsilons, flag


def spectral_weight(g, omega, epsilons, smearing=0.01):
    m = np.shape(epsilons)

    Omega = np.eye(m[0], dtype=np.complex) * (omega + smearing * 1j)
    G00 = LA.inv(Omega - epsilons)

    g00 = 0
    for i in range(g.num_wann):
        g00 += G00[i, i]

    surf_spectralN = -1.0 / np.pi * np.imag(g00)
    return surf_spectralN


def per_k(g, kpt):
    H00, H01 = construct_H00_H01(g, kpt)

    epsilon0 = H00
    alpha0 = H01
    beta0 = H01.conj().T
    surf_spectralN = []
    for i in range(g.dict['energy_div_num'] + 1):
        omega = g.eng_list[i]
        epsilon, epsilons, flag = iterative_method(alpha0, beta0, epsilon0, epsilon0, omega,
                                                   smearing=g.dict['smearing'], precision=g.dict['convergence'],
                                                   min_step=g.dict['minimum_iteration'],
                                                   max_step=g.dict['maximum_iteration'])

        surf_spectralN.append(spectral_weight(g, omega, epsilons, smearing=g.dict['smearing']))
        print(flag)

    return surf_spectralN
