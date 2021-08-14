# -*- coding: utf-8 -*-

from pathlib import Path

import pytest
import pydicom
import nibabel as nib
import numpy as np
from scipy.io import loadmat


# main paths
@pytest.fixture(scope='session')
def root_path() -> Path:
    return Path(__file__).parent


@pytest.fixture(scope='session')
def data_path(root_path) -> Path:
    return root_path / 'data'


# epi template
@pytest.fixture(scope='session')
def nii_image_1(data_path: Path) -> nib.nifti1.Nifti1Image:
    fp = data_path / 'fanon-0007-00006-000006-01.nii'
    return nib.load(fp, mmap=False)


# first test dcm
@pytest.fixture(scope='session')
def dcm_image(data_path: Path) -> np.array:
    fp = data_path / '001_000007_000006.dcm'
    return pydicom.dcmread(fp).pixel_array


# matlab resuls and settings
@pytest.fixture(scope='session')
def main_loop_data(data_path: Path) -> dict:
    fp = str(data_path / 'mainLoopData.mat')
    return loadmat(fp, squeeze_me=True)["mainLoopData"]


@pytest.fixture(scope='session')
def p_struct(data_path: Path) -> dict:
    fp = str(data_path / 'P.mat')
    return loadmat(fp, squeeze_me=True)["P"]


@pytest.fixture(scope='session')
def matlab_result(data_path: Path) -> np.array:
    fp = str(data_path / 'reslVol_matlab.mat')
    return loadmat(fp, squeeze_me=True)