import numpy as np
import pyshtools as pysh
from scipy.special import factorial
import warnings


def legendre(deg: int, x: np.ndarray) -> np.ndarray:
    """Computes Legendre functions of degree n.

    Notes:
        Equivalent to MATLAB legendre() with 'norm' normalization.
        We use the shtools '4pi' normalization, which comes out 2x
        larger than the MATLAB 'norm' scheme.
        The scipy.special.lpmv() function returns un-normalized values,
        hence using the pyshtools package.

    Args:
        deg (int): Degree of Legendre functions
        x (np.ndarray): Elements to evaluate
    """
    return (
        np.asarray([
            pysh.legendre.legendre_lm(deg, i, x, "4pi") for i in range(deg + 1)
        ])
        / 2
    )


def sph_to_farfield(Q, mmax, nmax, theta, phi):
    """
    This function computes the far fields E_theta and E_phi from the spherical modal coefficients Q,
    stored in single index format, suppressing the e{-jkr}/r term as is usual in FEKO.

    Parameters:
    Q (numpy array): Spherical modal coefficients
    mmax (int): Maximum azimuthal (phi) mode number
    nmax (int): Maximum elevation (theta) mode number
    theta (numpy array): 1D array of theta values in radians
    phi (numpy array): 1D array of phi values in radians

    Returns:
    E_th (numpy array): Far field in theta direction
    E_ph (numpy array): Far field in phi direction
    mode_counter (int): Mode counter
    """

    eta0 = 376.730313668  # post-2019 definition

    # Scale Q to align with FEKO usage
    # Also complex conjugate the coefficients to align with FEKO e^{j omega t} time convention
    Q = np.sqrt(8 * np.pi) * np.conj(Q)

    # Preallocate storage for speed

    num_modes = 2 * (nmax * (nmax + 1) + mmax - 1) + 2
    # print(f"{nmax} {mmax} {num_modes}")
    E_th_mode = np.zeros((len(phi), len(theta), num_modes + 1), dtype=np.complex128)
    E_ph_mode = np.zeros((len(phi), len(theta), num_modes + 1), dtype=np.complex128)

    mode_counter = 0
    for n in range(1, nmax + 1):
        # Calculate associated Legendre polynomials
        NP = legendre(n, np.cos(theta))

        # Add an extra row of zeros for the m+1 mode
        NP = np.pad(NP, ((0, 1), (0, 0)), mode="constant")

        for m in range(-n, n + 1):
            m_indx = -m  # The m-mode swaps due to the GRASP time convention
            if m >= 0:
                Sign = 1
            else:
                Sign = (-1) ** m

            CmnConstant = np.sqrt((n + abs(m) + 1) * (n - abs(m)))

            # Precompute m NP(x)/sin(x) and abs(m) NP(x)/sin(x) terms, including for special cases
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                NPdsm = NP[abs(m), :] / np.sin(theta) * m
                NPdsabsm = NP[abs(m), :] / np.sin(theta) * abs(m)

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

            # TE modes
            s = 1
            j = 2 * ((n + 1) * n + m_indx - 1) + s
            # print(f"TE j s n m {j} {s} {n} {m_indx}")
            for tt in range(len(theta)):
                phi_const = (np.exp(1j * m * phi) * Sign / np.sqrt(n * (n + 1))).T
                E_th_mode[:, tt, j - 1] = -(1j**n) * NPdsm[tt] * Q[j - 1]
                E_th_mode[:, tt, j - 1] = phi_const * E_th_mode[:, tt, j - 1]
                E_ph_mode[:, tt, j - 1] = -(1j ** (n + 1)) * (
                    NPdsabsm[tt] * Q[j - 1] * np.cos(theta[tt])
                    - CmnConstant * Q[j - 1] * NP[abs(m) + 1, tt]
                )
                E_ph_mode[:, tt, j - 1] = phi_const * E_ph_mode[:, tt, j - 1]
            mode_counter += 1

            # TM modes
            s = 2
            j = 2 * ((n + 1) * n + m_indx - 1) + s
            # print(f"TM j s n m {j} {s} {n} {m_indx}")
            for tt in range(len(theta)):
                phi_const = (np.exp(1j * m * phi) * Sign / np.sqrt(n * (n + 1))).T
                E_th_mode[:, tt, j - 1] = 1j**n * (
                    NPdsabsm[tt] * Q[j - 1] * np.cos(theta[tt])
                    - CmnConstant * Q[j - 1] * NP[abs(m) + 1, tt]
                )
                E_th_mode[:, tt, j - 1] = phi_const * E_th_mode[:, tt, j - 1]
                E_ph_mode[:, tt, j - 1] = 1j ** (n + 1) * NPdsm[tt] * Q[j - 1]
                E_ph_mode[:, tt, j - 1] = phi_const * E_ph_mode[:, tt, j - 1]
            mode_counter += 1

    # Now sum over all modes
    E_th = np.sum(E_th_mode, axis=2)
    E_ph = np.sum(E_ph_mode, axis=2)

    # Apply scaling needed by FEKO
    # Also apply 1/2 factor due to sphtools '4pi' vs. MATLAB 'norm' normalization
    sqrt_fac = np.sqrt(eta0 / (2 * np.pi))
    E_th = sqrt_fac * E_th
    E_ph = sqrt_fac * E_ph

    return E_th, E_ph, mode_counter