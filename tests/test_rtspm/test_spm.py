import numpy as np

from rtspm import spm_realign_rt
from rtspm import spm_reslice_rt


def img_2d_to_3d(img2d, xdim_img_number, ydim_img_number, dim3d):
    sl = 0
    vol3d = np.zeros(dim3d)
    for sy in range(0, ydim_img_number):
        for sx in range(0, xdim_img_number):
            if sl > dim3d[2]:
                break
            else:
                vol3d[:, :, sl] = img2d[sy * dim3d[0]: (sy + 1) * dim3d[0], sx * dim3d[0]: (sx + 1) * dim3d[0]]
            vol3d[:, :, sl] = np.rot90(vol3d[:, :, sl], 3)
            sl += 1

    return vol3d


def vol_3d_to_2d(vol3d, sl_nr_img2d_dimx, sl_nr_img2d_dimy, xdim_img_number, ydim_img_number, dim3d):
    sl = 0
    img_2d = np.zeros((ydim_img_number, xdim_img_number))

    for sy in range(0, sl_nr_img2d_dimy):
        for sx in range(0, sl_nr_img2d_dimx):
            if sl > dim3d[2]:
                break
            else:
                img_2d[sy * dim3d[1]:(sy + 1) * dim3d[1], sx * dim3d[0]:(sx + 1) * dim3d[0]] = np.rot90(vol3d[:, :, sl])
            sl += 1

    return img_2d


def get_mosaic_dim(dim3d):
    xdim_img_number = round(np.sqrt(dim3d[2]))
    tmp_dim = dim3d[2] - xdim_img_number ** 2

    if tmp_dim == 0:
        ydim_img_number = xdim_img_number
    else:
        if tmp_dim > 0:
            ydim_img_number = xdim_img_number
            xdim_img_number += 1
        else:
            xdim_img_number = xdim_img_number
            ydim_img_number = xdim_img_number

    img2d_dimx = xdim_img_number * dim3d[0]
    img2d_dimy = ydim_img_number * dim3d[0]

    return xdim_img_number, ydim_img_number, img2d_dimx, img2d_dimy


def test_spm(data_path, dcm_image, nii_image_1, p_struct, matlab_result):
    try:
        a0 = []
        x1 = []
        x2 = []
        x3 = []
        deg = []
        b = []
        r = [{'mat': np.array([]), 'dim': np.array([]), 'Vol': np.array([])} for _ in range(2)]
        dim_vol = np.array([74, 74, 36])

        r[0]["mat"] = nii_image_1.affine
        r[0]["dim"] = dim_vol.copy()
        tmp_vol = np.array(nii_image_1.get_fdata(), order='F')

        xdim_img_number, ydim_img_number, img2d_dimx, img2d_dimy = get_mosaic_dim(dim_vol)

        nr_zero_pad_vol = p_struct["nrZeroPadVol"].item()
        if p_struct["isZeroPadding"].item():
            r[0]["dim"][2] = r[0]["dim"][2] + nr_zero_pad_vol * 2
            r[0]["Vol"] = np.pad(tmp_vol, ((0, 0), (0, 0), (nr_zero_pad_vol, nr_zero_pad_vol)),
                                 'constant', constant_values=(0, 0))
        else:
            r[0]["Vol"] = tmp_vol

        ind_vol = 6
        dcm_data = np.array(dcm_image, dtype=float)

        r[1]["mat"] = nii_image_1.affine
        tmp_vol = img_2d_to_3d(dcm_data, xdim_img_number, ydim_img_number, dim_vol)

        if p_struct["isZeroPadding"].item():
            dim_vol[2] = dim_vol[2] + nr_zero_pad_vol * 2
            r[1]["Vol"] = np.pad(tmp_vol, ((0, 0), (0, 0), (nr_zero_pad_vol, nr_zero_pad_vol)), 'constant',
                                 constant_values=(0, 0))
        else:
            r[1]["Vol"] = tmp_vol

        r[1]["Vol"] = np.array(r[1]["Vol"], order='F')
        r[1]["dim"] = dim_vol

        flags_spm_realign = dict({'quality': .9, 'fwhm': 5, 'sep': 4, 'interp': 4, 'wrap': np.zeros((3, 1)),
                                  'rtm': 0, 'PW': '', 'lkp': np.array(range(0, 6))})
        flags_spm_reslice = dict({'quality': .9, 'fwhm': 5, 'sep': 4, 'interp': 4, 'wrap': np.zeros((3, 1)),
                                  'mask': 1, 'mean': 0, 'which': 2})

        nr_skip_vol = p_struct["nrSkipVol"].item()
        r, *_ = spm_realign_rt(r, flags_spm_realign, ind_vol, nr_skip_vol + 1,
                            a0, x1, x2, x3, deg, b)

        if p_struct["isZeroPadding"].item():
            tmp_resl_vol = spm_reslice_rt(r, flags_spm_reslice)
            resl_vol = tmp_resl_vol[:, :, nr_zero_pad_vol: -1 - nr_zero_pad_vol + 1]
            dim_vol[2] = dim_vol[2] - nr_zero_pad_vol * 2
        else:
            resl_vol = spm_reslice_rt(r, flags_spm_reslice)

        # resl_dic = {"reslVol_python": resl_vol}
        # savemat(data_path / "reslVol.mat", resl_dic)

        matlab_resl_vol = matlab_result["reslVol"]
        print('\n\nMSE = {:}\n'.format(((resl_vol - matlab_resl_vol) ** 2).mean()))

        assert True, "Done"
    except Exception as err:
        assert False, f"Error occurred: {repr(err)}"
