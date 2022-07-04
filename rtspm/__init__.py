# -*- coding: utf-8 -*-

from ._spm_inv_t_cdf import spm_inv_t_cdf
from ._spm_smoothing import spm_smooth
from ._spm_realign_rt import spm_realign_rt
from ._spm_reslice_rt import spm_reslice_rt
from ._spm_slice_vol import spm_slice_vol
from ._spm_matrix import spm_matrix
from ._spm_imatrix import spm_imatrix
from ._spm_setup import spm_setup
from .version import __version__


__all__ = [
    'spm_inv_t_cdf',
    'spm_smooth',
    'spm_realign_rt',
    'spm_reslice_rt',
    'spm_slice_vol',
    'spm_matrix',
    'spm_imatrix',
    'spm_setup',
    '__version__',
]
