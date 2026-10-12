"""ant_sph_tools.healpix -- Convert spherical harmonic coefficients between pyshtools and healpy.

Conventions (verified in ``tests/test_healpix.py``; see ``docs/healpy_vs_pyshtools.md``):

healpy
    Y_lm(theta, phi) = lambda_lm(cos theta) exp(i m phi),
    lambda_lm = sqrt((2l+1)/(4 pi) (l-m)!/(l+m)!) P_lm, *including* the
    Condon-Shortley phase (-1)^m.  Orthonormal: int |Y_lm|^2 dOmega = 1.
    Maps are real, so a_{l,-m} = (-1)^m conj(a_lm) and only m >= 0 is
    stored, in a 1D array ordered by m then l (``hp.Alm.getidx``).

pyshtools
    ``kind='real'``:  f = sum_{l, m>=0} [C_lm cos(m phi) + S_lm sin(m phi)] Pbar_lm
    with Pbar_lm carrying an extra sqrt(2 - delta_m0) for the 4pi, ortho and
    schmidt normalizations.  coeffs[0] = C_lm, coeffs[1] = S_lm.
    ``kind='complex'``: f = sum_{l, m} c_lm Pbar_l|m| exp(i m phi) with
    *no* sqrt(2) factor.  coeffs[0, l, m] = c_lm, coeffs[1, l, m] = c_{l,-m}.
    Default normalization is '4pi' and default csphase is +1 (Condon-Shortley
    phase *excluded*).

So pyshtools ``kind='complex', normalization='ortho', csphase=-1`` is exactly
the healpy basis.  Everything else needs a per-(l, m) rescaling, a (-1)^m
sign, and -- for real coefficients -- a_lm = (C_lm - i S_lm) / sqrt(2), m > 0.

healpy itself is not required: the alm index helpers below reproduce ``hp.Alm.getlm``,
``hp.Alm.getidx`` and ``hp.Alm.getsize`` for mmax == lmax.
"""

import numpy as np
import pyshtools as pysh
from scipy.special import gammaln

__all__ = [
    'alm_getidx',
    'alm_getlm',
    'alm_getlmax',
    'alm_getsize',
    'healpy_to_shtools',
    'norm_to_ortho_factor',
    'shtools_to_healpy',
]

_NORMS = ('4pi', 'ortho', 'schmidt', 'unnorm')


def alm_getsize(lmax: int) -> int:
    """Number of healpy alm coefficients for lmax (mmax = lmax). Same as ``hp.Alm.getsize``."""
    return (lmax + 1) * (lmax + 2) // 2


def alm_getlmax(size: int) -> int:
    """Infer lmax from an alm array length (mmax = lmax); -1 if invalid. Same as ``hp.Alm.getlmax``."""
    lmax = int(round((np.sqrt(8 * size + 1) - 3) / 2))
    return lmax if lmax >= 0 and alm_getsize(lmax) == size else -1


def alm_getidx(lmax: int, l, m):
    """Index of (l, m) in a healpy alm array (mmax = lmax). Same as ``hp.Alm.getidx``."""
    return m * (2 * lmax + 1 - m) // 2 + l


def alm_getlm(lmax: int):
    """Return (l, m) arrays for every entry of a healpy alm array. Same as ``hp.Alm.getlm``."""
    m = np.concatenate([np.full(lmax + 1 - mm, mm) for mm in range(lmax + 1)])
    l = np.concatenate([np.arange(mm, lmax + 1) for mm in range(lmax + 1)])
    return l, m


def norm_to_ortho_factor(l, m, normalization: str, kind: str) -> np.ndarray:
    """Factor F with coeff_ortho = F * coeff_<normalization> (same ``kind``).

    Args:
        l (array_like): Degree(s).
        m (array_like): Order(s); the sign of m is ignored.
        normalization: One of '4pi', 'ortho', 'schmidt', 'unnorm'.
        kind: 'real' or 'complex' (only matters for 'unnorm').

    Returns:
        np.ndarray of conversion factors, broadcast to shape of l and m.
    """
    l = np.asarray(l, dtype=float)
    m = np.abs(np.asarray(m, dtype=float))
    l, m = np.broadcast_arrays(l, m)
    if normalization == 'ortho':
        return np.ones_like(l)
    if normalization == '4pi':
        return np.full_like(l, np.sqrt(4 * np.pi))
    if normalization == 'schmidt':
        return np.sqrt(4 * np.pi / (2 * l + 1))
    if normalization == 'unnorm':
        f = np.sqrt(4 * np.pi / (2 * l + 1) * np.exp(gammaln(l + m + 1) - gammaln(l - m + 1)))
        if kind == 'real':
            f = f / np.sqrt(np.where(m == 0, 1.0, 2.0))
        return f
    raise ValueError(f'normalization must be one of {_NORMS}, got {normalization!r}')


def _cs_sign(m, csphase: int) -> np.ndarray:
    """Multiplier taking pyshtools coefficients to the healpy (CS-included) phase."""
    if csphase not in (1, -1):
        raise ValueError(f'csphase must be 1 or -1, got {csphase!r}')
    m = np.abs(np.asarray(m))
    return np.where((m % 2 == 1) & (csphase == 1), -1.0, 1.0)


def _unpack(clm, normalization, csphase, kind):
    """Return (coeffs, normalization, csphase, kind) from SHCoeffs or array."""
    if isinstance(clm, pysh.SHCoeffs):
        return clm.coeffs, clm.normalization, clm.csphase, clm.kind
    coeffs = np.asarray(clm)
    if coeffs.ndim != 3 or coeffs.shape[0] != 2 or coeffs.shape[1] != coeffs.shape[2]:
        raise ValueError(f'Expected shape (2, lmax+1, lmax+1), got {coeffs.shape}')
    if kind is None:
        kind = 'complex' if np.iscomplexobj(coeffs) else 'real'
    return coeffs, normalization or '4pi', 1 if csphase is None else csphase, kind


def shtools_to_healpy(
    clm,
    normalization: str = None,
    csphase: int = None,
    kind: str = None,
    lmax: int = None,
    complex_field: bool = False,
    rtol: float = 1e-8,
):
    """Convert pyshtools coefficients to a healpy ``alm`` array.

    Args:
        clm: ``pyshtools.SHCoeffs`` instance, or a ``(2, L+1, L+1)`` array as
            returned by ``pyshtools.expand`` routines / ``SHCoeffs.coeffs``.
        normalization: pyshtools normalization of an *array* input
            ('4pi' (pyshtools default), 'ortho', 'schmidt', 'unnorm').
            Ignored for SHCoeffs input (read from the object).
        csphase: pyshtools csphase of an array input (1 = excluded, default; -1 = included).
        kind: 'real' or 'complex' for an array input (default: inferred from dtype).
        lmax: Output lmax. Defaults to the input lmax; smaller truncates, larger zero-pads.
        complex_field: If True, return ``(alm_re, alm_im)`` -- healpy alms of the real
            and imaginary parts of the field -- so that
            ``f = hp.alm2map(alm_re) + 1j * hp.alm2map(alm_im)``. If False (default)
            the field must be real and a single alm array is returned.
        rtol: Tolerance on the imaginary part for the real-field check.

    Returns:
        np.ndarray (complex128) in healpy ordering, or a tuple of two if ``complex_field``.
    """
    coeffs, normalization, csphase, kind = _unpack(clm, normalization, csphase, kind)
    lmax_in = coeffs.shape[1] - 1
    lmax = lmax_in if lmax is None else int(lmax)

    l, m = alm_getlm(lmax)
    ok = l <= lmax_in  # entries beyond the input lmax stay zero
    li, mi = l[ok], m[ok]
    scale = norm_to_ortho_factor(li, mi, normalization, kind) * _cs_sign(mi, csphase)

    if kind == 'real':
        if np.iscomplexobj(coeffs) and np.any(coeffs.imag):
            raise ValueError("kind='real' but coefficients have non-zero imaginary part")
        C = coeffs[0].real[li, mi] * scale
        S = coeffs[1].real[li, mi] * scale
        a_pos = np.where(mi == 0, C, (C - 1j * S) / np.sqrt(2))
        alm = np.zeros(len(l), dtype=complex)
        alm[ok] = a_pos
        if complex_field:
            return alm, np.zeros_like(alm)
        return alm

    if kind != 'complex':
        raise ValueError(f"kind must be 'real' or 'complex', got {kind!r}")

    a_pos = coeffs[0][li, mi] * scale  # a_{l, m}
    a_neg = np.where(mi == 0, a_pos, coeffs[1][li, mi] * scale)  # a_{l, -m}
    sgn = (-1.0) ** mi
    # f* has coefficients (-1)^m conj(a_{l,-m}) in the healpy basis
    re = np.zeros(len(l), dtype=complex)
    im = np.zeros(len(l), dtype=complex)
    re[ok] = 0.5 * (a_pos + sgn * np.conj(a_neg))
    im[ok] = -0.5j * (a_pos - sgn * np.conj(a_neg))
    if complex_field:
        return re, im
    norm = max(np.max(np.abs(re)), np.finfo(float).tiny)
    if np.max(np.abs(im)) > rtol * norm:
        raise ValueError(
            'Coefficients describe a complex-valued field; healpy alms describe real maps. '
            'Use complex_field=True to get alms for the real and imaginary parts.'
        )
    return re


def healpy_to_shtools(
    alm,
    lmax: int = None,
    normalization: str = 'ortho',
    csphase: int = -1,
    kind: str = 'real',
    alm_imag=None,
) -> pysh.SHCoeffs:
    """Convert a healpy ``alm`` array to a ``pyshtools.SHCoeffs`` object.

    Args:
        alm: healpy alm array (real field, m >= 0 stored).
        lmax: lmax of ``alm`` (default: inferred assuming mmax == lmax).
        normalization (str): pyshtools normalization of the output ('ortho' (default), '4pi',
            'schmidt', 'unnorm').
        csphase (int): pyshtools Condon-Shortley phase of the output (-1 (default) or 1).
        kind (str): 'real' (default) or 'complex'.
        alm_imag: Optional alm of the imaginary part of a complex field
            (requires ``kind='complex'``).

    Returns:
        pyshtools.SHCoeffs
    """
    alm = np.asarray(alm, dtype=complex)
    lmax = alm_getlmax(len(alm)) if lmax is None else int(lmax)
    if lmax < 0 or alm_getsize(lmax) != len(alm):
        raise ValueError('alm size inconsistent with lmax (mmax < lmax not supported)')
    l, m = alm_getlm(lmax)
    inv = 1.0 / (norm_to_ortho_factor(l, m, normalization, kind) * _cs_sign(m, csphase))

    if kind == 'real':
        if alm_imag is not None:
            raise ValueError("alm_imag requires kind='complex'")
        coeffs = np.zeros((2, lmax + 1, lmax + 1))
        C = np.where(m == 0, alm.real, np.sqrt(2) * alm.real)
        S = np.where(m == 0, 0.0, -np.sqrt(2) * alm.imag)
        coeffs[0, l, m] = C * inv
        coeffs[1, l, m] = S * inv
    elif kind == 'complex':
        ai = np.zeros_like(alm) if alm_imag is None else np.asarray(alm_imag, dtype=complex)
        sgn = (-1.0) ** m
        pos = alm + 1j * ai  # c_{l, m}
        neg = sgn * np.conj(alm) + 1j * sgn * np.conj(ai)  # c_{l, -m}
        coeffs = np.zeros((2, lmax + 1, lmax + 1), dtype=complex)
        coeffs[0, l, m] = pos * inv
        coeffs[1, l, m] = np.where(m == 0, 0, neg * inv)
    else:
        raise ValueError(f"kind must be 'real' or 'complex', got {kind!r}")

    return pysh.SHCoeffs.from_array(coeffs, normalization=normalization, csphase=csphase, copy=False)
