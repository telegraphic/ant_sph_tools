"""legendre_norm.py - equivalent to MATLAB legendre() function with 'norm' normalization."""

import numpy as np
import pyshtools as pysh


def legendre_norm(deg: int, x: np.ndarray) -> np.ndarray:
    """Computes Legendre functions of degree n.

    Notes:
        Equivalent to MATLAB legendre() with 'norm' normalization.
        We use the SHTOOLS '4pi' normalization, which comes out 2x
        larger than the MATLAB 'norm' scheme.

        See https://shtools.github.io/SHTOOLS/complex-spherical-harmonics.html#4pi-normalized
        and https://au.mathworks.com/help/matlab/ref/legendre.html#f89-1002493

        SHTOOLS (and most of astro) uses (l, m) for degree and order, respectively.
        MATLAB (and engineering?) uses (n, m) for degree and order, respectively.

        The scipy.special.lpmv() function returns un-normalized values,
        hence using the pyshtools package.

    Args:
        deg (int): Degree of Legendre functions
        x (np.ndarray): Elements to evaluate

    Returns:
        f_leg = the associated Legendre functions of degree n and order m
    """
    f_leg = [pysh.legendre.legendre_lm(deg, i, x, '4pi') for i in range(deg + 1)]
    # Factor of 1/2 to convert 4pi normalization to MATLAB
    return 0.5 * np.asarray(f_leg)