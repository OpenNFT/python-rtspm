# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import betainc


def spm_inv_b_cdf(f, v, w, tol=1e-8):
    max_it = 1e4

    ad = np.array([2, 2, 2])
    rd = np.max(ad)
    az = np.array([[1, 1],
                   [1, 1],
                   [1, 1]])
    rs = np.max(az)
    xa = np.prod(az, 1) > 1

    x = np.zeros(rs)

    md = 0 <= f <= 1 and v > 0 and w > 0

    x[md and f == 1] = 1

    q = np.nonzero(md and 0 < f < 1)[0]
    if q.size == 0:
        return x

    if xa[0]:
        fq = f
        fq = fq.flatten()
    else:
        fq = f * np.ones((np.max(q.shape), 1))

    if xa[1]:
        vq = v
        vq = vq.flatten()
    else:
        vq = v * np.ones((np.max(q.shape), 1))

    if xa[2]:
        wq = w
        wq = wq.flatten()
    else:
        wq = w * np.ones((np.max(q.shape), 1))

    a = np.zeros((np.max(q.shape), 1))
    fa = a - fq
    b = np.ones((np.max(q.shape), 1))
    fb = b - fq
    i = 0
    xq = a + 0.5
    qq = np.array(range(0, (np.max(q.shape))))

    while qq.size > 0 and i < max_it:
        fxqq = betainc(vq[qq], wq[qq], xq[qq]) - fq[qq]
        mqq = fa[qq] * fxqq > 0

        a[qq[mqq[0]]] = xq[qq[mqq[0]]]
        fa[qq[mqq[0]]] = fxqq[mqq[0]]
        b[qq[not mqq]] = xq[qq[not mqq]]
        fb[qq[not mqq]] = fxqq[not mqq]
        xq[qq] = a[qq] + (b[qq] - a[qq]) / 2
        qq = qq[((b[qq] - a[qq]) > tol)[0]]

        i += 1

    x[q] = xq

    return x