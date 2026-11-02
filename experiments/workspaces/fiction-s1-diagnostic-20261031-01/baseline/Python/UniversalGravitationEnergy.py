import numpy as np

import Python_slovement.utils as utils

# Universal gravitation between each node and a set of fixed point-mass attractors,
# i.e. the potential energy p(x) = -G * m1 * m2 / ||x||_2 derived before:
#   val  : E = sum_i sum_j  -G * M_j * m_i / ||x_i - p_j||_2
#   grad : dE/dx_i = sum_j  G * M_j * m_i * (x_i - p_j) / ||x_i - p_j||^3
#   hess : d2E/dx_i^2 = sum_j  [ mu/r^3 * I - 3 * mu/r^5 * d d^T ],
#          with d = x_i - p_j, r = ||d||, mu = G * M_j * m_i
# The force F = -grad points toward each attractor with magnitude G * M_j * m_i / r^2.
# The exact Hessian is indefinite (radial eigenvalue -2*mu/r^3, tangential +mu/r^3),
# so project_hess_to_PSD flattens it for the projected Newton solver.

G = 6.674e-11                    # gravitational constant
attractors = []                  # [(mass, np.array([px, py])), ...], empty = module disabled
project_hess_to_PSD = True       # set False for the exact indefinite Hessian

def val(x, m):
    sum = 0.0
    for i in range(0, len(x)):
        for [M, p] in attractors:
            sum += -G * M * m[i] / np.linalg.norm(x[i] - p)
    return sum

def grad(x, m):
    g = np.array([[0.0, 0.0]] * len(x))
    for i in range(0, len(x)):
        for [M, p] in attractors:
            d = x[i] - p
            mu = G * M * m[i]
            g[i] += mu * d / np.linalg.norm(d) ** 3
    return g

def hess(x, m):
    IJV = [[], [], np.array([])]
    for i in range(0, len(x)):
        H_i = np.zeros((2, 2))
        for [M, p] in attractors:
            d = x[i] - p
            r = np.linalg.norm(d)
            mu = G * M * m[i]
            H_i += mu / r ** 3 * np.identity(2) - 3 * mu / r ** 5 * np.outer(d, d)
        if project_hess_to_PSD:
            H_i = utils.make_PSD(H_i)
        for rI in range(0, 2):
            for cI in range(0, 2):
                IJV[0].append(i * 2 + rI)
                IJV[1].append(i * 2 + cI)
                IJV[2] = np.append(IJV[2], H_i[rI, cI])
    return IJV
