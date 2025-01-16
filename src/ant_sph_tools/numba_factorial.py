"""numba_factorial.py - Provides a numba-compatible version of scipy.special.factorial."""
import numpy as np
from numba import njit
from scipy.special import factorial

# Create a lookup table by calling scipy
FACTORIAL_LUT = factorial(np.arange(0, 150))

@njit
def numba_factorial(n: int) -> float:
    """Numba-compatible factorial function.

    This uses a pre-computed LUT,

    Args:
        n (int): Number to compute factorial for.

    Returns:
        N (int): the factorial of n (N = n!)
    """
    return FACTORIAL_LUT[n]