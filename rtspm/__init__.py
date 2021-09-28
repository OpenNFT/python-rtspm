# -*- coding: utf-8 -*-

import logging

from rtspm.version import __version__
from ._spm_inv_t_cdf import spm_inv_t_cdf
from ._spm_smoothing import spm_smooth
from ._spm_realign import spm_realign
from ._spm_reslice import spm_reslice
from ._spm_imatrix import spm_imatrix

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = [
    'spm_inv_t_cdf',
    'spm_smooth',
    'spm_realign',
    'spm_reslice',
    'spm_imatrix',
    '__version__'
]
