# -*- coding: utf-8 -*-

import numpy as np

from ._spm_inv_b_cdf import spm_inv_b_cdf


def spm_inv_t_cdf(f, v):
    ad = np.array([2, 2])
    rd = np.max(ad)
    az = np.array([[1, 1],
                   [1, 1]])
    rs = np.max(az)
    xa = np.prod(az, 1) > 1

    x = np.zeros(rs)

    md = 0 <= f <= 1 and v > 0

    x[md and f == 0] = -np.inf
    x[md and f == 1] = +np.inf

    ml = (md and v == 1)
    if xa[0]:
        mlf = ml
    else:
        mlf = 1

    x[ml] = np.tan(np.pi * (f - 0.5))

    q = np.nonzero(md and f != 0.5 and 1 > f > 0 != v)[0]
    if q.size == 0:
        return x

    xqpos = f > 0.5
    bq = spm_inv_b_cdf(2 * (xqpos - (xqpos * 2 - 1) * f), v / 2, 0.5)
    x[q] = (xqpos * 2 - 1) * np.sqrt(v / bq - v)

    return x