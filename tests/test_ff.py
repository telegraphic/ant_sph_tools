"""Tests: sph_to_farfield."""

import numpy as np
import pylab as plt
from scipy.io import loadmat

from ant_sph_tools import read_sph, sph_to_farfield


def test_sph_to_farfield():
    """Test sph_to_farfield()."""
    sph_data = read_sph('sph/dipole_FarField1_299MHz.sph')
    mat_data = loadmat('tests/test_FFfromSph.mat')

    theta = np.deg2rad(np.arange(0, 91))
    phi = np.deg2rad(np.arange(0, 361))

    E_th, E_ph, mode_counter = sph_to_farfield(
        sph_data['Q'], sph_data['mmax'], sph_data['nmax'], theta, phi
    )

    assert np.allclose(E_ph, mat_data['E_ph'])
    assert np.allclose(E_th, mat_data['E_th'])

    plt.figure('E_theta', figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.imshow(np.abs(E_th), aspect='auto')
    plt.colorbar()

    plt.subplot(1, 2, 2)
    plt.imshow(np.abs(mat_data['E_th']), aspect='auto')
    plt.colorbar()

    plt.figure('E_phi', figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.imshow(np.abs(E_ph), aspect='auto')
    plt.colorbar()

    plt.subplot(1, 2, 2)
    plt.imshow(np.abs(mat_data['E_ph']), aspect='auto')
    plt.colorbar()

    plt.show()


if __name__ == '__main__':
    test_sph_to_farfield()
