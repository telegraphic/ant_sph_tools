"""ant_sph_tools.pysh -- Spherical mode (.sph) data as pyshtools coefficients.

The far field of a spherical wave expansion (Hansen 1988; TICRA/GRASP .sph files) is a
tangential vector field on the sphere. Any such field can be written with two scalar
potentials (a Helmholtz-Hodge decomposition):

    E(theta, phi) = grad_s U  -  r_hat x grad_s V

i.e.

    E_theta = dU/dtheta + (1/sin theta) dV/dphi
    E_phi   = (1/sin theta) dU/dphi - dV/dtheta

The TM modes (s=2) map onto the poloidal potential U, and the TE modes (s=1) onto the
toroidal potential V. Expanding U and V in ordinary (scalar) spherical harmonics turns
the Q_smn coefficients into two sets of pyshtools coefficients:

    u_nm = sqrt(8 pi eta0) i^n     conj(Q_{2,-m,n}) / sqrt(n (n+1))
    v_nm = sqrt(8 pi eta0) i^(n+1) conj(Q_{1,-m,n}) / sqrt(n (n+1))

with u_nm, v_nm in pyshtools complex, 'ortho', csphase=1 convention. The sqrt(8 pi),
conjugation, m -> -m swap and sqrt(eta0 / 2 pi) all follow FFfromSph.m (Davidson & Sutinjo
2024), so the far fields match ``sph_to_farfield`` (and the MATLAB code) exactly.

The far field is then evaluated from U and V using pyshtools' normalized associated Legendre
functions, with pole-safe recurrences for dP/dtheta and m P / sin(theta), vectorised over m
and theta.
"""

import numpy as np
import pyshtools as pysh

from .io import read_sph

ETA0 = 376.730313668  # post-2019 definition, as in FFfromSph.m

__all__ = [
    'q_to_shcoeffs',
    'read_sph_shcoeffs',
    'shcoeffs_to_q',
    'shcoeffs_to_farfield',
    'sph_to_farfield_pysh',
]


def _mode_table(mmax: int, nmax: int):
    """Return arrays (n, m, j_te) of all modes with 1 <= n <= nmax, |m| <= min(n, mmax).

    j_te is the zero-based index of the TE (s=1) coefficient in the compressed Q array for
    mode (n, m); the TM (s=2) coefficient is at j_te + 1. Following FFfromSph.m, mode m is
    read from compressed index -m, due to the GRASP e^{-i omega t} time convention.
    """
    n = np.concatenate([np.full(2 * min(nn, mmax) + 1, nn) for nn in range(1, nmax + 1)])
    m = np.concatenate([np.arange(-min(nn, mmax), min(nn, mmax) + 1) for nn in range(1, nmax + 1)])
    j_te = 2 * (n * (n + 1) - m - 1)
    return n, m, j_te


def _scale(n: np.ndarray) -> np.ndarray:
    """Common scale factor sqrt(8 pi eta0) / sqrt(n (n+1))."""
    return np.sqrt(8 * np.pi * ETA0) / np.sqrt(n * (n + 1))


def q_to_shcoeffs(Q: np.ndarray, mmax: int, nmax: int, normalization: str = 'ortho', csphase: int = -1):
    """Convert spherical mode coefficients Q into pyshtools coefficients of U and V.

    Args:
        Q (np.ndarray): Spherical mode coefficients in compressed index order (from ``read_sph``).
        mmax (int): Maximum azimuthal mode number.
        nmax (int): Maximum polar mode number.
        normalization (str): pyshtools normalization of the output ('ortho', '4pi', 'schmidt',
            'unnorm'). Default 'ortho'.
        csphase (int): pyshtools Condon-Shortley phase convention of the output. Default -1
            (included), which together with 'ortho' is the healpy convention.

    Returns:
        U (pysh.SHCoeffs): Complex coefficients of the poloidal potential (TM modes, s=2).
        V (pysh.SHCoeffs): Complex coefficients of the toroidal potential (TE modes, s=1).
    """
    Q = np.asarray(Q)
    n, m, j_te = _mode_table(mmax, nmax)
    if j_te.max() + 1 >= len(Q):
        raise ValueError(f'Q has {len(Q)} entries; too few for mmax={mmax}, nmax={nmax}')

    k = _scale(n)
    u = k * (1j**n) * np.conj(Q[j_te + 1])
    v = k * (1j ** (n + 1)) * np.conj(Q[j_te])

    cu = np.zeros((2, nmax + 1, nmax + 1), dtype=complex)
    cv = np.zeros_like(cu)
    i = (m < 0).astype(int)  # pyshtools stores m < 0 in coeffs[1, l, |m|]
    cu[i, n, np.abs(m)] = u
    cv[i, n, np.abs(m)] = v

    U = pysh.SHCoeffs.from_array(cu, normalization='ortho', csphase=1)
    V = pysh.SHCoeffs.from_array(cv, normalization='ortho', csphase=1)
    if (normalization, csphase) != ('ortho', 1):
        U = U.convert(normalization=normalization, csphase=csphase)
        V = V.convert(normalization=normalization, csphase=csphase)
    return U, V


def read_sph_shcoeffs(fn: str, normalization: str = 'ortho', csphase: int = -1) -> dict:
    """Read a TICRA/GRASP .sph file and return pyshtools coefficients.

    Args:
        fn (str): File name of the .sph file.
        normalization (str): pyshtools normalization of the output. Default 'ortho'.
        csphase (int): pyshtools csphase of the output. Default -1 (healpy-compatible).

    Returns:
        sph_data (dict): Dictionary with entries:
            U (pysh.SHCoeffs): Poloidal potential coefficients (TM modes).
            V (pysh.SHCoeffs): Toroidal potential coefficients (TE modes).
            mmax (int): Maximum azimuthal mode number.
            nmax (int): Maximum polar mode number.
    """
    d = read_sph(fn)
    U, V = q_to_shcoeffs(d['Q'], d['mmax'], d['nmax'], normalization=normalization, csphase=csphase)
    return {'U': U, 'V': V, 'mmax': d['mmax'], 'nmax': d['nmax']}


def _to_ortho_cs1(clm: pysh.SHCoeffs) -> np.ndarray:
    """Return the (2, L+1, L+1) complex coefficient array in 'ortho', csphase=1 convention."""
    if clm.kind == 'real':
        clm = clm.convert(kind='complex')
    if (clm.normalization, clm.csphase) != ('ortho', 1):
        clm = clm.convert(normalization='ortho', csphase=1)
    return np.asarray(clm.coeffs, dtype=complex)


def shcoeffs_to_q(U: pysh.SHCoeffs, V: pysh.SHCoeffs, mmax: int = None, nmax: int = None) -> np.ndarray:
    """Inverse of ``q_to_shcoeffs``: convert U, V back into compressed Q coefficients.

    Args:
        U (pysh.SHCoeffs): Poloidal potential coefficients (any pyshtools convention).
        V (pysh.SHCoeffs): Toroidal potential coefficients (any pyshtools convention).
        mmax (int): Maximum azimuthal mode number (default: lmax).
        nmax (int): Maximum polar mode number (default: lmax).

    Returns:
        Q (np.ndarray): Spherical mode coefficients in compressed index order, as from ``read_sph``.
    """
    cu, cv = _to_ortho_cs1(U), _to_ortho_cs1(V)
    lmax = max(cu.shape[1], cv.shape[1]) - 1
    nmax = lmax if nmax is None else nmax
    mmax = nmax if mmax is None else mmax
    cu = _pad(cu, nmax)
    cv = _pad(cv, nmax)

    n, m, j_te = _mode_table(mmax, nmax)
    i = (m < 0).astype(int)
    k = _scale(n)
    Q = np.zeros(2 * (nmax * (nmax + 1) + mmax - 1) + 2, dtype=complex)
    Q[j_te + 1] = np.conj(cu[i, n, np.abs(m)] / (k * 1j**n))
    Q[j_te] = np.conj(cv[i, n, np.abs(m)] / (k * 1j ** (n + 1)))
    return Q


def _pad(c: np.ndarray, lmax: int) -> np.ndarray:
    """Truncate or zero-pad a (2, L+1, L+1) coefficient array to lmax."""
    out = np.zeros((2, lmax + 1, lmax + 1), dtype=c.dtype)
    L = min(lmax, c.shape[1] - 1)
    out[:, : L + 1, : L + 1] = c[:, : L + 1, : L + 1]
    return out


def _legendre_terms(lmax: int, theta: np.ndarray):
    """Compute dPbar_lm/dtheta and m Pbar_lm / sin(theta) for 0 <= m <= l <= lmax.

    Pbar_lm are the pyshtools 'ortho', complex-normalized (cnorm=1), csphase=1 Legendre
    functions. Both quantities are evaluated with recurrences that avoid dividing by
    sin(theta), so they are exact at the poles.

    Returns:
        A (np.ndarray): dPbar_lm/dtheta, shape (Ntheta, lmax+1, lmax+1).
        B (np.ndarray): m Pbar_lm / sin(theta), shape (Ntheta, lmax+1, lmax+1).
    """
    z = np.cos(np.asarray(theta, dtype=float))
    P = np.zeros((len(z), lmax + 2, lmax + 2))  # extra row/col of zeros for m+1 lookups
    for t, zz in enumerate(z):
        P[t, : lmax + 1, : lmax + 1] = pysh.legendre.legendre(
            lmax, zz, normalization='ortho', csphase=1, cnorm=1
        )

    l = np.arange(lmax + 1)[:, None].astype(float)
    m = np.arange(lmax + 1)[None, :].astype(float)
    Pm_minus = np.zeros((len(z), lmax + 1, lmax + 1))
    Pm_minus[:, :, 1:] = P[:, : lmax + 1, :lmax]  # P_{l, m-1}
    Pm_plus = P[:, : lmax + 1, 1 : lmax + 2]  # P_{l, m+1}

    # dP_lm/dtheta = 1/2 [sqrt((l+m)(l-m+1)) P_{l,m-1} - sqrt((l-m)(l+m+1)) P_{l,m+1}],  m >= 1
    # dP_l0/dtheta = -sqrt(l(l+1)) P_{l,1}
    a = np.sqrt(np.clip((l + m) * (l - m + 1), 0, None))
    b = np.sqrt(np.clip((l - m) * (l + m + 1), 0, None))
    A = 0.5 * (a * Pm_minus - b * Pm_plus)
    A[:, :, 0] = -np.sqrt(l[:, 0] * (l[:, 0] + 1)) * P[:, : lmax + 1, 1]

    # m P_lm / sin(theta) = 1/2 sqrt((2l+1)/(2l-1))
    #                       [sqrt((l+m)(l+m-1)) P_{l-1,m-1} + sqrt((l-m)(l-m-1)) P_{l-1,m+1}],  m >= 1
    Pl1_minus = np.zeros((len(z), lmax + 1, lmax + 1))
    Pl1_plus = np.zeros((len(z), lmax + 1, lmax + 1))
    Pl1_minus[:, 1:, 1:] = P[:, :lmax, :lmax]  # P_{l-1, m-1}
    Pl1_plus[:, 1:, :] = P[:, :lmax, 1 : lmax + 2]  # P_{l-1, m+1}
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.where(l > 0, 0.5 * np.sqrt((2 * l + 1) / (2 * l - 1)), 0.0)
    c1 = np.sqrt(np.clip((l + m) * (l + m - 1), 0, None))
    c2 = np.sqrt(np.clip((l - m) * (l - m - 1), 0, None)) * (m <= l - 2)
    B = c * (c1 * Pl1_minus + c2 * Pl1_plus)
    B[:, :, 0] = 0.0
    B *= m <= l
    A *= m <= l
    return A, B


def shcoeffs_to_farfield(U: pysh.SHCoeffs, V: pysh.SHCoeffs, theta: np.ndarray, phi: np.ndarray):
    """Compute the far fields E_theta, E_phi from pyshtools potential coefficients U and V.

    Args:
        U (pysh.SHCoeffs): Poloidal potential coefficients (any pyshtools convention).
        V (pysh.SHCoeffs): Toroidal potential coefficients (any pyshtools convention).
        theta (np.ndarray): 1D array of theta (colatitude) values in radians.
        phi (np.ndarray): 1D array of phi values in radians.

    Returns:
        E_th (np.ndarray): Far field in theta direction, shape (len(phi), len(theta)).
        E_ph (np.ndarray): Far field in phi direction, shape (len(phi), len(theta)).
    """
    theta = np.atleast_1d(np.asarray(theta, dtype=float))
    phi = np.atleast_1d(np.asarray(phi, dtype=float))
    cu, cv = _to_ortho_cs1(U), _to_ortho_cs1(V)
    lmax = max(cu.shape[1], cv.shape[1]) - 1
    cu, cv = _pad(cu, lmax), _pad(cv, lmax)

    A, B = _legendre_terms(lmax, theta)

    # Sum over degree l for each order: F[t, |m|]
    def lsum(c, X):
        return np.einsum('lm,tlm->tm', c, X)

    up, un, vp, vn = cu[0], cu[1], cv[0], cv[1]
    k = np.arange(lmax + 1)
    sgn = (-1.0) ** k
    # m >= 0:  E_th <- u A + i v B,  E_ph <- i u B - v A
    Fth_p = lsum(up, A) + 1j * lsum(vp, B)
    Fph_p = 1j * lsum(up, B) - lsum(vp, A)
    # m = -k < 0 (Pbar_{l,-k} = (-1)^k Pbar_lk):  E_th <- (-1)^k (u A - i v B),  E_ph <- (-1)^k (-i u B - v A)
    Fth_n = sgn * (lsum(un, A) - 1j * lsum(vn, B))
    Fph_n = sgn * (-1j * lsum(un, B) - lsum(vn, A))

    eip = np.exp(1j * np.outer(phi, k))  # (Nphi, lmax+1)
    E_th = eip @ Fth_p.T + np.conj(eip[:, 1:]) @ Fth_n[:, 1:].T
    E_ph = eip @ Fph_p.T + np.conj(eip[:, 1:]) @ Fph_n[:, 1:].T
    return E_th, E_ph


def sph_to_farfield_pysh(Q: np.ndarray, mmax: int, nmax: int, theta: np.ndarray, phi: np.ndarray):
    """Drop-in replacement for ``sph_to_farfield`` that goes via pyshtools coefficients.

    Args:
        Q (np.ndarray): Spherical modal coefficients, in compressed index format.
        mmax (int): Maximum azimuthal (phi) mode number.
        nmax (int): Maximum elevation (theta) mode number.
        theta (np.ndarray): 1D array of theta values in radians.
        phi (np.ndarray): 1D array of phi values in radians.

    Returns:
        E_th (np.ndarray): Far field in theta direction, shape (len(phi), len(theta)).
        E_ph (np.ndarray): Far field in phi direction, shape (len(phi), len(theta)).
        mode_counter (int): Number of (TE + TM) modes used.
    """
    U, V = q_to_shcoeffs(Q, mmax, nmax, normalization='ortho', csphase=1)
    E_th, E_ph = shcoeffs_to_farfield(U, V, theta, phi)
    mode_counter = 2 * len(_mode_table(mmax, nmax)[0])
    return E_th, E_ph, mode_counter
