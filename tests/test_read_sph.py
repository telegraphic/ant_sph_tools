"""Tests: test file read."""

import numpy as np
from scipy.io import loadmat

from ant_sph_tools.io import read_sph


def test_read_sph():
    """Test read_sph function.

    Compares against stored MATLAB output.
    """
    sph_data = read_sph('sph/dipole_FarField1_299MHz.sph')

    # Load MATLAB workspace
    mat_data = loadmat('tests/test_ReadSphModes.mat')

    assert np.allclose(sph_data['Q'], mat_data['Q'])
    assert np.allclose(sph_data['j'], mat_data['j'])
    assert sph_data['mmax'] == mat_data['mmax']
    assert sph_data['nmax'] == mat_data['nmax']


if __name__ == '__main__':
    test_read_sph()
