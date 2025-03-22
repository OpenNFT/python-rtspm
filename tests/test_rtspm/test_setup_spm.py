import numpy as np
from scipy.io import savemat

from rtspm import _spm_setup


def test_setup_spm(p_cont_struct, protocols_inds_cont_struct, spm_cont_struct,
                   p_int_struct, protocols_inds_int_struct, spm_int_struct):

    tr = p_cont_struct["TR"].item()
    mean_vol_templ = p_cont_struct["meanVolTemplate"].item()
    nscan = p_cont_struct["NrOfVolumes"].item() - p_cont_struct["nrSkipVol"].item()

    offsets = [protocols_inds_cont_struct["offsets"]]
    first_inds = [protocols_inds_cont_struct["first_nf_inds"]]
    prot_names = [protocols_inds_cont_struct["prot_names"]]

    spm_python = _spm_setup.spm_setup(tr, nscan, mean_vol_templ, offsets, first_inds, prot_names)

    spm_matlab_x = spm_cont_struct["x"]
    spm_matlab_x0 = spm_cont_struct["x0"]

    spm_python_x = spm_python["xX_x"]
    spm_python_x0 = spm_python["xX_K"]["x0"]

    savemat("spm_python_cont.mat", {"x_p": spm_python_x, "x0_p": spm_python_x0})

    np.testing.assert_almost_equal(spm_matlab_x, spm_python_x, decimal=7, err_msg="Not equal")
    np.testing.assert_almost_equal(spm_matlab_x0, spm_python_x0, decimal=7, err_msg="Not equal")

    tr = p_int_struct["TR"].item()
    mean_vol_templ = p_int_struct["meanVolTemplate"].item()
    nscan = p_int_struct["NrOfVolumes"].item() - p_int_struct["nrSkipVol"].item()

    offsets = list(protocols_inds_int_struct["offsets"])
    first_inds = protocols_inds_int_struct["first_nf_inds"]
    prot_names = protocols_inds_int_struct["prot_names"]

    spm_python = _spm_setup.spm_setup(tr, nscan, mean_vol_templ, offsets, first_inds, prot_names)

    spm_matlab_x = spm_int_struct["x"]
    spm_matlab_x0 = spm_int_struct["x0"]

    spm_python_x = spm_python["xX_x"]
    spm_python_x0 = spm_python["xX_K"]["x0"]

    savemat("spm_python_int.mat", {"x_p": spm_python_x, "x0_p": spm_python_x0})

    np.testing.assert_almost_equal(spm_matlab_x, spm_python_x, decimal=7, err_msg="Not equal")
    np.testing.assert_almost_equal(spm_matlab_x0, spm_python_x0, decimal=7, err_msg="Not equal")







