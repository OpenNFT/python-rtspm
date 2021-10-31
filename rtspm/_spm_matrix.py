# -*- coding: utf-8 -*-

import numpy as np


def spm_matrix(p):
    if p.size == 3:
        a = np.eye(4)
        a[0:3, 3] = p[:]
        return a

    q = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0])
    p = np.append(p, q[p.size:12])

    t = np.array([[1, 0, 0, p[0]],
                  [0, 1, 0, p[1]],
                  [0, 0, 1, p[2]],
                  [0, 0, 0, 1]])

    r1 = np.array([[1, 0, 0, 0],
                   [0, np.cos(p[3]), np.sin(p[3]), 0],
                   [0, -np.sin(p[3]), np.cos(p[3]), 0],
                   [0, 0, 0, 1]])

    r2 = np.array([[np.cos(p[4]), 0, np.sin(p[4]), 0],
                   [0, 1, 0, 0],
                   [-np.sin(p[4]), 0, np.cos(p[4]), 0],
                   [0, 0, 0, 1]])

    r3 = np.array([[np.cos(p[5]), np.sin(p[5]), 0, 0],
                   [-np.sin(p[5]), np.cos(p[5]), 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]])

    r = r1 @ r2 @ r3

    z = np.array([[p[6], 0, 0, 0],
                  [0, p[7], 0, 0],
                  [0, 0, p[8], 0],
                  [0, 0, 0, 1]])

    s = np.array([[1, p[9], p[10], 0],
                  [0, 1, p[11], 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])

    a = t @ r @ z @ s

    return a