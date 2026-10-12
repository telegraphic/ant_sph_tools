## ant_sph_tools

This repository implements the methods in  Davidson, D. B., & Sutinjo, A. T. (2024). "Efficient storage of embedded element patterns for low frequency radio telescopes". _Radio Science_, 59, e2024RS008080.
DOI: [10.1029/2024RS008080](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024RS008080)

Specifically, this is a Python implementation of the MATLAB code provided in the accompanying dataset "MATLAB files for reconstructing antenna patterns from spherical mode data" DOI: [10.25917/n05j-th75](https://doi.org/10.25917/N05J-TH75)

The code provides two main functions: `read_sph()`, for reading spherical mode data in TIRCA/GRASP .sph format, and `sph_to_farfield()`, which computes the far-field pattern from the .sph data:

```python
def read_sph(fn: str) -> dict:
    """Reads in spherical mode data in TIRCA .sph format.

    Args:
        fn (str): File name of the .sph file

    Returns:
        sph_data (dict): Dictionary with following entries:
            Q (numpy array): Spherical mode coefficients
            j (numpy array): Compressed mode indices
            j3D (numpy array): 3D mode indices
            mmax (int): Maximum azimuthal mode number
            nmax (int): Maximum elevation mode number
    """

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
```

### pyshtools coefficients

`ant_sph_tools.pysh` re-expresses the spherical mode data as two sets of [pyshtools](https://shtools.github.io/SHTOOLS/) coefficients. The far field is a tangential vector field, so it can be written with two scalar potentials (a Helmholtz decomposition):

```
E = grad U - r_hat x grad V
E_theta = dU/dtheta + (1/sin theta) dV/dphi
E_phi   = (1/sin theta) dU/dphi - dV/dtheta
```

The TM modes (s=2) give the poloidal potential `U` and the TE modes (s=1) the toroidal potential `V`:

```
u_nm = sqrt(8 pi eta0) i^n     conj(Q_{2,-m,n}) / sqrt(n(n+1))
v_nm = sqrt(8 pi eta0) i^(n+1) conj(Q_{1,-m,n}) / sqrt(n(n+1))
```

(pyshtools complex, `'ortho'`, `csphase=1`; the scaling, conjugation and m -> -m swap follow `FFfromSph.m`.)

```python
from ant_sph_tools import read_sph_shcoeffs, shcoeffs_to_farfield, sph_to_farfield_pysh

d = read_sph_shcoeffs('tests/sph/dipole_FarField1_299MHz.sph')   # {'U', 'V', 'mmax', 'nmax'}
E_th, E_ph = shcoeffs_to_farfield(d['U'], d['V'], theta, phi)    # shape (len(phi), len(theta))

# drop-in replacement for sph_to_farfield (same arguments and outputs)
E_th, E_ph, mode_counter = sph_to_farfield_pysh(Q, mmax, nmax, theta, phi)
```

`U` and `V` are complex `pysh.SHCoeffs`; the default convention (`'ortho'`, `csphase=-1`) is the healpy convention, and any other pyshtools normalization/phase can be requested. `shcoeffs_to_q()` converts back to `Q`. The far-field evaluation uses pyshtools Legendre functions with pole-safe recurrences, is vectorised, and agrees with the MATLAB test data to ~1e-13 (about 60x faster than `sph_to_farfield` for nmax = 99). Tests: `tests/test_pysh.py`.

### healpy <-> pyshtools

`ant_sph_tools.healpix` converts spherical harmonic coefficients between pyshtools (any normalization, Condon-Shortley phase, real or complex kind) and [healpy](https://healpy.readthedocs.io/) `alm` arrays. healpy is only needed to *use* the alm, not to convert them.

```python
from ant_sph_tools import healpy_to_shtools, shtools_to_healpy

alm = shtools_to_healpy(clm)                                       # SHCoeffs or (2, L+1, L+1) array
alm_re, alm_im = shtools_to_healpy(d['U'], complex_field=True)     # complex fields -> two alm arrays
clm = healpy_to_shtools(alm, normalization='4pi', csphase=1, kind='real')
```

The two libraries differ in phase convention (healpy includes the Condon-Shortley phase, pyshtools excludes it by default), normalization (orthonormal vs 4pi), real vs complex harmonics (real pyshtools coefficients map to `(C - iS)/sqrt(2)` for m > 0, not `C + iS`), and storage. See [docs/healpy_vs_pyshtools.md](docs/healpy_vs_pyshtools.md) for the full comparison and conversion formulas. Tests: `tests/test_healpix.py`.

### Requirements

The code requires numpy, scipy and the [pyshtools](https://shtools.github.io/SHTOOLS/) package. The healpy conversion tests also require healpy (skipped if it is not installed).

### Tests

The code output has been checked for equivalence against its MATLAB counterpart; test data exported from MATLAB is provided in the `tests` directory.

### Development with pixi

```
pixi install
```

Tasks:

```
pixi run test       # run tests
pixi run cov        # run code coverage
pixi run lint       # run ruff linter
pixi run jupyter    # run jupyterlab
```
