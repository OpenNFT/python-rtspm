# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import gammaln


def spm_setup(tr, nscan, mean_vol_template, offsets, first_inds, prot_names):
    spm = {
        "xY_RT": tr / 1000,
        "nscan": nscan,
        "xBF_T": 16,
        "xBF_T0": 1,
        "xBF_UNITS": 'scans',
        "xBF_Volterra": 1,
        "xBF_name": 'hrf',
        "xBF_length": 32,
        "xBF_order": 1,
        "xBF_dt": tr / 16000,
        "xX_K_HParam": 128,
        "sess_C_C": None,
        "sess_C_name": None,
        "xM_TH": np.ones((nscan, 1)) * mean_vol_template
    }

    fmri_t = spm["xBF_T"]
    fmri_t0 = spm["xBF_T0"]

    dt = spm["xBF_dt"]
    bf, p = spm_hrf(dt, fmri_t)

    spm["xBF_length"] = bf.shape[0] * dt
    spm["xBF_order"] = bf.shape[1]

    spm["xBF_bf"] = bf
    n = spm["nscan"]

    u = {"name": [], "ons": [], "dur": [], "u": []}
    for i in range(len(prot_names)):
        u["name"].append(prot_names[i])
        u["ons"].append(first_inds[i])
        u["dur"].append((np.diff(offsets[i].T, axis=0).T + 1).squeeze())

    u = spm_get_ons(u, n, fmri_t, dt, spm["xBF_T"] * spm["xBF_dt"])

    x, fc = spm_volterra(u, spm["xBF_bf"])

    if x.size > 0:
        x = x[np.array(range(0, n)) * fmri_t + fmri_t0 + 32, :]

    spm["sess_row"] = np.array(range(0, n), ndmin=2)
    spm["xX_x"] = np.hstack((x, np.ones((n, 1))))

    k = {"h_param": spm["xX_K_HParam"], "row": spm["sess_row"], "RT": spm["xY_RT"]}

    spm["xX_K"] = spm_filter(k)

    return spm


def spm_filter(k):
    k_1 = k["row"].shape[1]
    n = np.fix(2 * (k_1 * k["RT"]) / k["h_param"] + 1).astype(np.int32)
    x0 = spm_dctmtx(k_1, n)
    k["x0"] = x0[:, 1:]

    return k


def spm_dctmtx(n, k):
    n_1 = np.array(range(0, n)).T
    c = np.zeros((n, k))
    c[:, 0] = (np.ones((n, 1)) / np.sqrt(n)).squeeze()
    for j in range(1, k):
        c[:, j] = np.sqrt(2 / n) * np.cos(np.pi * (2 * n_1 + 1) * j / (2 * n))

    return c


def spm_volterra(u, bf):
    x = np.array([])
    fc = {"i": [], "name": [], "p": []}
    for i in range(len(u["name"])):
        ind = []
        ip = []
        for k in range(u["u"][i].shape[1]):
            for p in range(bf.shape[1]):
                x_1 = u["u"][i][:, k]
                d = np.array(range(0, len(x_1)))
                x_1 = np.convolve(x_1, bf[:, p])
                x_1 = x_1[d]
                if x.any():
                    x_1 = np.array(x_1, ndmin=2).T
                    x = np.hstack((x, x_1))
                else:
                    x = np.array(x_1, ndmin=2).T

                ind.append(x.shape[1])
                ip.append(k)

        fc["i"].append(ind)
        fc["name"].append(u["name"][i])
        fc["p"].append(ip)

    return x, fc


def spm_get_ons(u, k, fmri_t, dt, tr):
    for i in range(len(u["name"])):

        # input ons is Python indexes, starting from zero
        ons = u["ons"][i] + 1
        dur = u["dur"][i]

        uu = np.array(ons ** 0, ndmin=2).T
        ton = np.array(np.round(ons * tr / dt) + 33, ndmin=2, dtype=np.int32)
        tof = np.array(np.round(dur * tr / dt) + ton + 1, ndmin=2, dtype=np.int32)
        sf = np.zeros((k * fmri_t + 128, uu.shape[1]))
        ton = np.max(ton, axis=0)
        tof = np.max(tof, axis=0)
        for j in range(0, len(ton)):
            if sf.shape[0] > ton[j]:
                sf[ton[j], :] = sf[ton[j], :] + uu[j, :]
            if sf.shape[0] > tof[j]:
                sf[tof[j], :] = sf[tof[j], :] - uu[j, :]
        sf = np.cumsum(sf)
        sf = np.array(sf[0:(k * fmri_t + 32)], ndmin=2).T

        u["u"].append(sf)

    return u


def spm_hrf(rt, fmri_t):
    p = np.array([6, 16, 1, 1, 6, 0, 32])

    dt = rt / fmri_t
    u = np.array(range(0, (p[6] / dt).__ceil__() + 1), ndmin=2)
    hrf = spm_Gpdf(u, p[0] / p[2], dt / p[2]) - spm_Gpdf(u, p[1] / p[3], dt / p[3]) / p[4]
    hrf = hrf[np.array(range(0, (p[6] / rt).__floor__() + 1), ndmin=2) * fmri_t]
    hrf = hrf.T / np.sum(hrf, axis=None)

    return hrf, p


def spm_Gpdf(x, h, l):
    # Probability Density Function (PDF) of Gamma distribution

    ac = [list(x.shape), [1, 1], [1, 1]]

    rc = np.max(ac)
    xa = np.array(np.prod(ac, axis=1) > 1, dtype=np.int32)

    f = np.zeros((rc,))

    md = (np.ones(x.shape).any() and h > 0 and l > 0)
    ml = np.where(md and x.any() and h < 1)
    f[ml] = np.inf
    ml = np.where(md and x.any() and h == 1)
    f[ml] = l

    q = np.where(md and x > 0)[1]
    if q.size == 0:
        return
    if xa[0]:
        qx = q
    else:
        qx = 0

    f[q] = np.exp(
        ((h - 1) * np.log(x[0, qx]) + h * np.log(l) - l * x[0, qx] - gammaln(h))
    )

    return f
