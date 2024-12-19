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
```
