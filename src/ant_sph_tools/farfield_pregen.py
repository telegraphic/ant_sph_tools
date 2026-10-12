"""farfield.py -- Convert sph file to farfield E_theta E_phi."""

import warnings
import time
import h5py
import numpy as np
from tqdm import tqdm

from .legendre_norm import legendre_norm as legendre
from .numba_factorial import numba_factorial as factorial


def _sph_to_farfield(Q, mmax, nmax, theta, phi):
    """sph_to_farfield compute kernel. Called by sph_to_farfield()."""
    # Scale Q to align with FEKO usage
    # Also complex conjugate the coefficients to align with FEKO e^{j omega t} time convention
    Q = np.sqrt(8 * np.pi) * np.conj(Q)

    # Preallocate storage for speed
    num_modes = 2 * (nmax * (nmax + 1) + mmax - 1) + 2

    with h5py.File('modes.h5', 'r') as h5:
        E_th_modes = h5['E_theta_modes'][:]
        E_ph_modes = h5['E_phi_modes'][:]

    mode_counter = 0

    js = []
    for n in tqdm(range(1, nmax + 1)):
        for m in range(-n, n + 1):
            # Compressed mode number j: uses -m (sign swap) due to the GRASP time convention
            j = 2 * ((n + 1) * n - m - 1)
            js.append(j)
            js.append(j+1)

    j = np.array(js)
    print(j)

    t0 = time.time()
    E_th = np.sum(E_th_modes[:, :, j] * Q[j], axis=2)
    E_ph = np.sum(E_ph_modes[:, :, j] * Q[j], axis=2)
    t1 = time.time()
    print(f"{t1-t0:.2f}s mult")


    # Scaling needed by FEKO
    # Also apply 1/2 factor due to sphtools '4pi' vs. MATLAB 'norm' normalization
    eta0 = 376.730313668  # post-2019 definition
    sqrt_fac = np.sqrt(eta0 / (2 * np.pi))

    # Now sum over all modes
    E_th = sqrt_fac * E_th
    E_ph = sqrt_fac * E_ph

    return E_th, E_ph, mode_counter


def regenerate(Q, mmax, nmax, theta, phi):
    """This function computes the far fields E_theta and E_phi from spherical modal coefficients Q.

    Args:
        Q (numpy array): Spherical modal coefficients, in compressed index format.
        mmax (int): Maximum azimuthal (phi) mode number
        nmax (int): Maximum elevation (theta) mode number
        theta (numpy array): 1D array of theta values in radians
        phi (numpy array): 1D array of phi values in radians

    Returns:
        E_th (numpy array): Far field in theta direction
        E_ph (numpy array): Far field in phi direction
        mode_counter (int): Mode counter
    """
    E_th, E_ph, mode_counter = _sph_to_farfield(Q, mmax, nmax, theta, phi)

    return E_th, E_ph, mode_counter
