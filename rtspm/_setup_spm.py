import numpy as np
from scipy.special import gammaln


def setup_spm(TR, nscan, mean_vol_template):
    THR = 0.5

    SPM = {"xY_RT": TR / 1000,
           "nscan": nscan,
           "xBF_T": 16,
           "xBF_T0": 1,
           "xBF_UNITS": 'scans',
           "xBF_Volterra": 1,
           "xBF_name": 'hrf',
           "xBF_length": 32,
           "xBF_order": 1,
           "xBF_dt": TR / 16000,
           "xX_K_HParam": 128,
           "sess_C_C": None,
           "sess_C_name": None,
           "xM_TH": np.ones((nscan, 1)) * mean_vol_template
           }

    # spm_fmri_ui

    # spm_fMRI_design

    fmri_t = SPM["xBF_T"]
    fmri_t0 = SPM["xBF_T0"]

    # spm_get_bf

    dt = SPM["xBF_dt"]
    bf, p = spm_hrf(dt, fmri_t)







def spm_hrf(dt, fmri_t):

    p = np.array([6, 16, 1, 1, 6, 0, 32])

    u = np.array(range(0, (p[6] / dt).__ceil__() + 1), ndmin=2)
    hrf = spm_Gpdf(u,p[0]/p[2],dt/p[2]) - spm_Gpdf(u,p[1]/p[3],dt/p[3]/p[4])
    hrf = hrf[ np.array(range(0, (p[6] / dt).__floor__() + 1), ndmin=2)*fmri_t ]
    hrf = hrf.T / np.sum(hrf,axis=None)

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
        ((h - 1) * np.log(x[0,qx]) + h * np.log(l) - l * x[0,qx] - gammaln(h))
    )

    return f
