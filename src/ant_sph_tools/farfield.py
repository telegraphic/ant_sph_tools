"""farfield.py -- Convert sph file to farfield E_theta E_phi."""

import warnings

import numpy as np

from .legendre_norm import legendre_norm as legendre
from .numba_factorial import numba_factorial as factorial


def m_mode_sign(m: int) -> int:
    """Return sign (+/- 1) required in calcs for m-modes."""
    if m >= 0:
        return 1
    else:
        return (-1) ** m


def fix_special_cases(n, m, theta, NPdsm, NPdsabsm):
    """Fix calculations for special cases."""
    indx0 = np.where(np.isclose(theta, 0))[0]
    if len(indx0) > 0:  # special case for theta=0
        if abs(m) != 1:
            NPdsm[indx0] = 0
            NPdsabsm[indx0] = 0
        else:
            Cmn = np.sqrt(
                (2 * n + 1) / 2 * factorial(n - abs(m)) / factorial(n + abs(m))
            )
            NPdsm[indx0] = Cmn * np.sign(m) * n * (n + 1) / 2
            NPdsabsm[indx0] = Cmn * n * (n + 1) / 2

    indx1 = np.where(np.isclose(theta, np.pi))[0]
    if len(indx1) > 0:  # special case for theta=180 deg
        if abs(m) != 1:
            NPdsm[indx1] = 0
            NPdsabsm[indx1] = 0
        else:
            Cmn = np.sqrt(
                (2 * n + 1) / 2 * factorial(n - abs(m)) / factorial(n + abs(m))
            )
            NPdsm[indx1] = (-1) ** (n + 1) * Cmn * np.sign(m) * n * (n + 1) / 2
            NPdsabsm[indx1] = (-1) ** (n + 1) * Cmn * n * (n + 1) / 2

    return NPdsm, NPdsabsm


def _sph_to_farfield(Q, mmax, nmax, theta, phi, legendre_list):
    """sph_to_farfield compute kernel. Called by sph_to_farfield()."""
    # Scale Q to align with FEKO usage
    # Also complex conjugate the coefficients to align with FEKO e^{j omega t} time convention
    Q = np.sqrt(8 * np.pi) * np.conj(Q)

    # Preallocate storage for speed
    num_modes = 2 * (nmax * (nmax + 1) + mmax - 1) + 2
    # print(f"{nmax} {mmax} {num_modes}")
    E_th_mode = np.zeros((len(phi), len(theta)), dtype=np.complex64)
    E_ph_mode = np.zeros((len(phi), len(theta)), dtype=np.complex64)

    # precompute reused values
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)

    mode_counter = 0
    for ii, n in enumerate(range(1, nmax + 1)):
        # Load precomputed Legendre polynomials
        NP = legendre_list[ii]

        # Precompute NP(x)/sin(x)
        NP_sin_theta = NP / sin_theta

        # Precompute some values that depend on m
        m_arr = np.arange(-n, n+1)
        CmnConstant_NP_arr = (np.sqrt((n + abs(m_arr) + 1) * (n - abs(m_arr))))[:, np.newaxis] * NP[abs(m_arr) + 1, :]

        NPdsm_arr = NP_sin_theta[abs(m_arr), :] * m_arr[:, np.newaxis]
        NPdsabsm_arr = NP_sin_theta[abs(m_arr), :] * abs(m_arr[:, np.newaxis])

        for cc, m in enumerate(range(-n, n + 1)):
            # Compute m NP(x)/sin(x) and abs(m) NP(x)/sin(x) terms, including for special cases
            NPdsm = NPdsm_arr[cc]
            NPdsabsm = NPdsabsm_arr[cc]
            NPdsm, NPdsabsm = fix_special_cases(n, m, theta, NPdsm, NPdsabsm)

            # Compressed mode number j: uses -m (sign swap) due to the GRASP time convention
            j = 2 * ((n + 1) * n - m - 1)

            # Compute sums that are reused
            phi_const = (np.exp(1j * m * phi) * m_mode_sign(m) / np.sqrt(n * (n + 1)))
            phi_const = phi_const[:, np.newaxis]

            Cmn_cos_factor = NPdsabsm * cos_theta - CmnConstant_NP_arr[cc]
            phi_const_cmn = phi_const * Cmn_cos_factor
            phi_const_NP = phi_const * NPdsm

            # TE modes
            E_th_mode += -(1j**n) * Q[j] * phi_const_NP
            E_ph_mode += -(1j ** (n + 1)) * Q[j] * phi_const_cmn
            mode_counter += 1

            # TM modes
            j += 1
            E_th_mode += 1j**n * Q[j] * phi_const_cmn
            E_ph_mode += 1j ** (n + 1) * Q[j] * phi_const_NP
            mode_counter += 1

    # Scaling needed by FEKO
    # Also apply 1/2 factor due to sphtools '4pi' vs. MATLAB 'norm' normalization
    eta0 = 376.730313668  # post-2019 definition
    sqrt_fac = np.sqrt(eta0 / (2 * np.pi))

    # Now sum over all modes
    E_th = sqrt_fac * E_th_mode
    E_ph = sqrt_fac * E_ph_mode

    return E_th, E_ph, mode_counter


def _calculate_legendre(nmax: int, theta: np.ndarray) -> list:
    """Compute required legendre polynomials."""
    legendre_list = []
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    for n in range(1, nmax + 1):
        # Calculate associated Legendre polynomials
        NP = legendre(n, cos_theta)

        # Add an extra row of zeros for the m+1 mode
        NP = np.pad(NP, ((0, 1), (0, 0)), mode='constant')
        legendre_list.append(NP)

    # Uncomment if using numba
    #legendre_list = nb.typed.List(legendre_list)
    return legendre_list


def sph_to_farfield(Q, mmax, nmax, theta, phi):
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
    # Precompute legendre polynomials
    legendre_list = _calculate_legendre(nmax, theta)

    # Call the main numba kernel
    with warnings.catch_warnings():
        # Suppress divide by zeros warnings
        warnings.simplefilter('ignore')
        E_th, E_ph, mode_counter = _sph_to_farfield(Q, mmax, nmax, theta, phi, legendre_list)

    return E_th, E_ph, mode_counter
