# -*- coding: utf-8 -*-

import numpy as np

from ._spm_matrix import spm_matrix


def spm_imatrix(m):
    r = m[0:3, 0:3]
    c = np.linalg.cholesky(r.T @ r)

    p = np.append(m[0:3, 3].T, np.zeros(3, ))
    p = np.append(p, np.diag(c).T)
    p = np.append(p, np.zeros(3, ))

    if np.linalg.det(r) < 0:
        p[6] = -p[6]

    c = np.linalg.solve(np.diag(np.diag(c)), c)
    p[9:12] = c.flatten()[[3, 6, 7]]
    r0 = spm_matrix(np.append(np.zeros(6, ), p[6:12]))
    r0 = r0[0:3, 0:3]
    r1 = np.linalg.solve(r0.T, r.T).T

    def rang(x):
        return np.minimum(np.maximum(x, -1), 1)

    p[4] = np.arcsin(rang(r1[0, 2]))
    if (np.abs(p[4]) - np.pi / 2) ** 2 < 1e-9:
        p[3] = 0
        p[5] = np.arctan2(-rang(r1[1, 0]), rang(-r1[2, 0] / r1[0, 2]))
    else:
        c = np.cos(p[4])
        p[3] = np.arctan2(rang(r1[1, 2] / c), rang(r1[2, 2] / c))
        p[5] = np.arctan2(rang(r1[0, 1] / c), rang(r1[0, 0] / c))

    return p