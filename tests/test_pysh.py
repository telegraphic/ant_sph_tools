"""Tests: pyshtools coefficients for .sph data, and the pyshtools far-field evaluation."""

import os

import numpy as np
import pytest
from scipy.io import loadmat

from ant_sph_tools import (
    read_sph,
    read_sph_shcoeffs,
    shcoeffs_to_farfield,
    shcoeffs_to_q,
    sph_to_farfield,
    sph_to_farfield_pysh,
)
from ant_sph_tools.pysh import q_to_shcoeffs

MATLAB_TESTS = [
    'dipole_FarField1_299MHz',
    'hertzian_dipole_FarField1_299MHz',
    'hertzian_x_dipole_FarField1_299MHz',
    'hertzian_y_dipole_FarField1_299MHz',
    'hertzian_xy_dipole_FarField1_299MHz',
    'hertzian_x_dip_array_FarField2_299MHz',
    'hertzian_z_dip_array_FarField1_299MHz',
]
# Committed reference file with m = 0, 2, 4 and both TE and TM modes
REF_SPH = 'tests/sph/hertzian_z_dip_array_FarField1_299MHz.sph'
# Larger AAVS element (nmax = 20); optional, test is skipped if the file is not present
AAVS_SPH = 'tests/sph-aavs/AAVS1d5_single_element_rev1_110MHz_FarField1(0).sph'

THETA = np.deg2rad(np.arange(0, 91))
PHI = np.deg2rad(np.arange(0, 361))


def assert_field_close(a, b, rtol=1e-10):
    """Assert two fields agree to rtol relative to the peak amplitude."""
    np.testing.assert_allclose(a, b, rtol=0, atol=rtol * np.abs(b).max())


@pytest.mark.parametrize('bn', MATLAB_TESTS)
def test_matlab_equivalence(bn):
    """sph_to_farfield_pysh reproduces FFfromSph.m output, including at theta = 0."""
    sph_data = read_sph(f'tests/sph/{bn}.sph')
    mat_data = loadmat(f'tests/matlab/test_{bn}.mat')
    E_th, E_ph, _ = sph_to_farfield_pysh(sph_data['Q'], sph_data['mmax'], sph_data['nmax'], THETA, PHI)
    assert E_th.shape == mat_data['E_th'].shape
    assert_field_close(E_th, mat_data['E_th'])
    assert_field_close(E_ph, mat_data['E_ph'])


@pytest.mark.parametrize('fn', [REF_SPH, AAVS_SPH])
def test_matches_sph_to_farfield(fn):
    """Agrees with the loop implementation, including a realistic AAVS element (nmax = 20) if present."""
    if not os.path.exists(fn):
        pytest.skip(f'{fn} not available')
    sph_data = read_sph(fn)
    theta = np.deg2rad(np.arange(0, 181))
    args = (sph_data['Q'], sph_data['mmax'], sph_data['nmax'], theta, PHI)
    E_th0, E_ph0, mc0 = sph_to_farfield(*args)
    E_th, E_ph, mc = sph_to_farfield_pysh(*args)
    assert mc == mc0
    # sph_to_farfield accumulates in complex64, so only ~1e-6 agreement is expected
    assert_field_close(E_th, E_th0, rtol=1e-5)
    assert_field_close(E_ph, E_ph0, rtol=1e-5)


@pytest.mark.parametrize('normalization', ['ortho', '4pi', 'schmidt'])
@pytest.mark.parametrize('csphase', [1, -1])
def test_read_sph_shcoeffs_any_convention(normalization, csphase):
    """Far field is independent of the pyshtools convention chosen for U and V."""
    d = read_sph_shcoeffs(REF_SPH, normalization=normalization, csphase=csphase)
    assert d['U'].normalization == normalization and d['U'].csphase == csphase
    assert d['U'].kind == 'complex' and d['U'].lmax == d['nmax']
    E_th, E_ph = shcoeffs_to_farfield(d['U'], d['V'], THETA, PHI)
    sph_data = read_sph(REF_SPH)
    E_th0, E_ph0, _ = sph_to_farfield_pysh(sph_data['Q'], sph_data['mmax'], sph_data['nmax'], THETA, PHI)
    assert_field_close(E_th, E_th0)
    assert_field_close(E_ph, E_ph0)


@pytest.mark.parametrize('fn', [f'tests/sph/{bn}.sph' for bn in MATLAB_TESTS])
def test_round_trip_q(fn):
    """Q -> (U, V) -> Q is lossless."""
    sph_data = read_sph(fn)
    U, V = q_to_shcoeffs(sph_data['Q'], sph_data['mmax'], sph_data['nmax'], normalization='4pi', csphase=1)
    Q = shcoeffs_to_q(U, V, mmax=sph_data['mmax'], nmax=sph_data['nmax'])
    np.testing.assert_allclose(Q, sph_data['Q'], rtol=1e-12, atol=1e-12 * np.abs(sph_data['Q']).max())


def test_helmholtz_decomposition():
    """E_theta = dU/dtheta + dV/dphi / sin(theta), E_phi = dU/dphi / sin(theta) - dV/dtheta.

    U and V are evaluated independently with pyshtools SHCoeffs.expand, and differentiated
    numerically, away from the poles.
    """
    d = read_sph_shcoeffs(REF_SPH)
    theta = np.deg2rad(np.linspace(5, 175, 35))
    phi = np.deg2rad(np.linspace(0, 355, 72))
    th, ph = np.meshgrid(theta, phi)
    h = 1e-5

    def f(clm, t, p):
        return clm.expand(lat=90 - np.rad2deg(t), lon=np.rad2deg(p), degrees=True)

    def d_theta(clm):
        return (f(clm, th + h, ph) - f(clm, th - h, ph)) / (2 * h)

    def d_phi(clm):
        return (f(clm, th, ph + h) - f(clm, th, ph - h)) / (2 * h)

    U, V = d['U'], d['V']
    E_th_fd = d_theta(U) + d_phi(V) / np.sin(th)
    E_ph_fd = d_phi(U) / np.sin(th) - d_theta(V)
    E_th, E_ph = shcoeffs_to_farfield(U, V, theta, phi)
    assert_field_close(E_th, E_th_fd, rtol=1e-6)
    assert_field_close(E_ph, E_ph_fd, rtol=1e-6)


def test_hertzian_z_dipole_is_pure_tm_n1_m0():
    """A z-directed Hertzian dipole is the single TM mode n=1, m=0: E_phi = 0, E_theta ~ sin(theta)."""
    d = read_sph_shcoeffs('tests/sph/hertzian_dipole_FarField1_299MHz.sph')
    U, V = d['U'].coeffs, d['V'].coeffs
    peak = np.abs(U).max()
    assert np.abs(U[0, 1, 0]) == pytest.approx(peak)
    U_other = U.copy()
    U_other[0, 1, 0] = 0
    assert np.abs(U_other).max() < 1e-12 * peak
    assert np.abs(V).max() < 1e-12 * peak

    theta = np.deg2rad(np.arange(0, 181))
    E_th, E_ph = shcoeffs_to_farfield(d['U'], d['V'], theta, PHI)
    assert np.abs(E_ph).max() < 1e-12 * np.abs(E_th).max()
    pattern = np.abs(E_th) / np.abs(E_th).max()
    np.testing.assert_allclose(pattern, np.broadcast_to(np.sin(theta), pattern.shape), atol=1e-12)


@pytest.mark.parametrize('bn', MATLAB_TESTS)
def test_pole_continuity(bn):
    """At theta = 0 and 180 deg the Cartesian field (E_x, E_y) must not depend on phi."""
    d = read_sph_shcoeffs(f'tests/sph/{bn}.sph')
    E_th, E_ph = shcoeffs_to_farfield(d['U'], d['V'], np.array([0.0, np.pi]), PHI)
    c, s = np.cos(PHI)[:, None], np.sin(PHI)[:, None]
    cos_th = np.array([1.0, -1.0])
    E_x = E_th * c * cos_th - E_ph * s
    E_y = E_th * s * cos_th + E_ph * c
    scale = max(np.abs(E_th).max(), np.abs(E_ph).max(), 1e-300)
    for E in (E_x, E_y):
        assert np.abs(E - E[0]).max() < 1e-10 * scale


def test_mmax_less_than_nmax():
    """Truncating mmax drops |m| > mmax modes only."""
    sph_data = read_sph(REF_SPH)
    U, V = q_to_shcoeffs(sph_data['Q'], 2, sph_data['nmax'])
    m = np.arange(U.lmax + 1)
    assert np.all(U.coeffs[:, :, m > 2] == 0) and np.all(V.coeffs[:, :, m > 2] == 0)
    assert np.any(U.coeffs[:, :, 2] != 0) and np.any(V.coeffs[:, :, 2] != 0)


def test_q_too_short_raises():
    """Asking for more modes than Q holds raises a clear error."""
    sph_data = read_sph('tests/sph/hertzian_dipole_FarField1_299MHz.sph')
    with pytest.raises(ValueError, match='too few'):
        q_to_shcoeffs(sph_data['Q'], sph_data['mmax'] + 2, sph_data['nmax'] + 2)
