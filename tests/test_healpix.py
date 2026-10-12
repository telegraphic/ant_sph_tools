"""Tests: pyshtools <-> healpy spherical harmonic coefficient conversion."""

import numpy as np
import pyshtools as pysh
import pytest

from ant_sph_tools import (
    healpy_to_shtools,
    read_sph_shcoeffs,
    shcoeffs_to_farfield,
    shtools_to_healpy,
)
from ant_sph_tools.healpix import alm_getidx, alm_getlm, alm_getlmax, alm_getsize

hp = pytest.importorskip('healpy')

LMAX = 12
NSIDE = 16
NORMS = ['4pi', 'ortho', 'schmidt', 'unnorm']
CSPHASES = [1, -1]


@pytest.fixture(scope='module')
def pix():
    """Colatitude/longitude (rad) and lat/lon (deg) of every HEALPix pixel centre."""
    theta, phi = hp.pix2ang(NSIDE, np.arange(hp.nside2npix(NSIDE)))
    return theta, phi, 90 - np.degrees(theta), np.degrees(phi)


def random_real_clm(lmax, normalization, csphase, seed=0):
    """Random real-kind SHCoeffs describing a real field, in the requested convention."""
    rng = np.random.default_rng(seed)
    c = rng.normal(size=(2, lmax + 1, lmax + 1))
    l = np.arange(lmax + 1)[:, None]
    m = np.arange(lmax + 1)[None, :]
    c[:, m > l] = 0
    c[1, :, 0] = 0
    clm = pysh.SHCoeffs.from_array(c, normalization='ortho', csphase=-1)
    return clm.convert(normalization=normalization, csphase=csphase, kind='real')


def random_complex_clm(lmax, normalization, csphase, seed=1):
    """Random complex-kind SHCoeffs describing a *complex-valued* field."""
    rng = np.random.default_rng(seed)
    c = rng.normal(size=(2, lmax + 1, lmax + 1)) + 1j * rng.normal(size=(2, lmax + 1, lmax + 1))
    l = np.arange(lmax + 1)[:, None]
    m = np.arange(lmax + 1)[None, :]
    c[:, m > l] = 0
    c[1, :, 0] = 0
    return pysh.SHCoeffs.from_array(c, normalization=normalization, csphase=csphase)


def eval_shtools(clm, lat, lon):
    """Evaluate pyshtools coefficients at (lat, lon) in degrees."""
    return clm.expand(lat=lat, lon=lon, degrees=True)


# ----------------------------------------------------------------------------------
# 1. Pin down the underlying conventions (these are the facts the note relies on)
# ----------------------------------------------------------------------------------
def test_healpy_includes_condon_shortley_phase(pix):
    """Setting healpy a_11 = 1 gives 2 Re(Y_11) with the (-1)^m phase included."""
    theta, phi, _, _ = pix
    alm = np.zeros(hp.Alm.getsize(1), complex)
    alm[hp.Alm.getidx(1, 1, 1)] = 1.0
    m = hp.alm2map(alm, NSIDE, lmax=1)
    Y11 = -np.sqrt(3 / (8 * np.pi)) * np.sin(theta) * np.exp(1j * phi)  # with (-1)^m
    np.testing.assert_allclose(m, 2 * Y11.real, atol=1e-12)  # a_{1,-1} term gives 2 Re


@pytest.mark.parametrize('m_', [0, 1, 2, 3])
def test_shtools_complex_ortho_csm1_is_healpy_basis(pix, m_):
    """SHTOOLS complex/ortho/csphase=-1 harmonics equal healpy Y_lm for +m and -m."""
    _, _, lat, lon = pix
    l_ = 3
    for sign in (+1, -1):
        c = np.zeros((2, l_ + 1, l_ + 1), complex)
        c[0 if sign > 0 else 1, l_, m_] = 1
        if sign < 0 and m_ == 0:
            continue
        ysh = pysh.expand.MakeGridPointC(c, lat, lon, norm=4, csphase=-1)
        # healpy basis function via alm2map of a complex field: Y = map(a) + i map(b)
        re, im = shtools_to_healpy(c, normalization='ortho', csphase=-1, complex_field=True)
        yhp = hp.alm2map(re, NSIDE, lmax=l_) + 1j * hp.alm2map(im, NSIDE, lmax=l_)
        np.testing.assert_allclose(ysh, yhp, atol=1e-12)
        # ... and the converted coefficients are just a single unit entry (no rescaling)
        if sign > 0 and m_ > 0:
            expected = np.zeros(hp.Alm.getsize(l_), complex)
            expected[hp.Alm.getidx(l_, l_, m_)] = 0.5
            np.testing.assert_allclose(re, expected, atol=1e-15)


def test_real_harmonics_carry_sqrt2():
    """Real ortho C_11 = 1 is sqrt(2) * Re(healpy Y_11) (up to CS sign), not Re(Y_11)."""
    c = np.zeros((2, 2, 2))
    c[0, 1, 1] = 1
    alm = shtools_to_healpy(c, normalization='ortho', csphase=-1, kind='real')
    assert alm[hp.Alm.getidx(1, 1, 1)] == pytest.approx(1 / np.sqrt(2))


# ----------------------------------------------------------------------------------
# 2. Converted coefficients reproduce the same map, for every convention
# ----------------------------------------------------------------------------------
@pytest.mark.parametrize('normalization', NORMS)
@pytest.mark.parametrize('csphase', CSPHASES)
@pytest.mark.parametrize('kind', ['real', 'complex'])
def test_real_field_maps_agree(pix, normalization, csphase, kind):
    """Converted alm reproduce the pyshtools field for every convention."""
    _, _, lat, lon = pix
    clm = random_real_clm(LMAX, normalization, csphase)
    if kind == 'complex':
        clm = clm.convert(kind='complex')
    f_sh = eval_shtools(clm, lat, lon)
    alm = shtools_to_healpy(clm)
    f_hp = hp.alm2map(alm, NSIDE, lmax=LMAX)
    np.testing.assert_allclose(np.real(f_sh), f_hp, atol=1e-10 * np.abs(f_hp).max())
    if kind == 'complex':
        assert np.abs(np.imag(f_sh)).max() < 1e-10 * np.abs(f_hp).max()


@pytest.mark.parametrize('normalization', NORMS)
@pytest.mark.parametrize('csphase', CSPHASES)
def test_array_input_matches_shcoeffs_input(normalization, csphase):
    """Raw arrays plus convention kwargs give the same result as SHCoeffs."""
    clm = random_real_clm(LMAX, normalization, csphase)
    a1 = shtools_to_healpy(clm)
    a2 = shtools_to_healpy(clm.coeffs, normalization=normalization, csphase=csphase)
    np.testing.assert_allclose(a1, a2)


@pytest.mark.parametrize('normalization', NORMS)
@pytest.mark.parametrize('csphase', CSPHASES)
@pytest.mark.parametrize('kind', ['real', 'complex'])
def test_independent_of_input_convention(normalization, csphase, kind):
    """Any pyshtools convention of the same field gives the same healpy alm."""
    ref = random_real_clm(LMAX, 'ortho', -1).convert(kind='complex')
    expected = shtools_to_healpy(ref)
    # complex/ortho/csphase=-1 is the healpy basis, so the alm are literally coeffs[0][l, m]
    l, m = hp.Alm.getlm(LMAX)
    np.testing.assert_allclose(expected, ref.coeffs[0][l, m], atol=1e-14)
    other = ref.convert(normalization=normalization, csphase=csphase, kind=kind)
    np.testing.assert_allclose(shtools_to_healpy(other), expected, atol=1e-10)


@pytest.mark.parametrize('normalization', NORMS)
@pytest.mark.parametrize('csphase', CSPHASES)
def test_complex_field(pix, normalization, csphase):
    """Complex fields convert to a pair of alm (real and imaginary parts)."""
    _, _, lat, lon = pix
    clm = random_complex_clm(LMAX, normalization, csphase)
    f_sh = eval_shtools(clm, lat, lon)
    re, im = shtools_to_healpy(clm, complex_field=True)
    f_hp = hp.alm2map(re, NSIDE, lmax=LMAX) + 1j * hp.alm2map(im, NSIDE, lmax=LMAX)
    np.testing.assert_allclose(f_sh, f_hp, atol=1e-10 * np.abs(f_hp).max())
    with pytest.raises(ValueError, match='complex-valued field'):
        shtools_to_healpy(clm)


# ----------------------------------------------------------------------------------
# 3. Round trips and lmax handling
# ----------------------------------------------------------------------------------
@pytest.mark.parametrize('normalization', NORMS)
@pytest.mark.parametrize('csphase', CSPHASES)
@pytest.mark.parametrize('kind', ['real', 'complex'])
def test_round_trip(pix, normalization, csphase, kind):
    """Alm -> pyshtools -> alm is lossless and preserves the map."""
    _, _, lat, lon = pix
    rng = np.random.default_rng(3)
    alm = rng.normal(size=hp.Alm.getsize(LMAX)) + 1j * rng.normal(size=hp.Alm.getsize(LMAX))
    l, m = hp.Alm.getlm(LMAX)
    alm[m == 0] = alm[m == 0].real
    clm = healpy_to_shtools(alm, normalization=normalization, csphase=csphase, kind=kind)
    assert clm.normalization == normalization and clm.csphase == csphase and clm.kind == kind
    np.testing.assert_allclose(shtools_to_healpy(clm), alm, atol=1e-10)
    f_hp = hp.alm2map(alm, NSIDE, lmax=LMAX)
    np.testing.assert_allclose(np.real(eval_shtools(clm, lat, lon)), f_hp, atol=1e-9)


def test_round_trip_complex_field(pix):
    """Complex-field alm pair -> pyshtools recovers the original coefficients."""
    _, _, lat, lon = pix
    clm = random_complex_clm(LMAX, '4pi', 1)
    re, im = shtools_to_healpy(clm, complex_field=True)
    back = healpy_to_shtools(re, alm_imag=im, normalization='4pi', csphase=1, kind='complex')
    np.testing.assert_allclose(back.coeffs, clm.coeffs, atol=1e-12)


def test_lmax_truncate_and_pad():
    """The lmax argument truncates or zero-pads."""
    clm = random_real_clm(LMAX, '4pi', 1)
    full = shtools_to_healpy(clm)
    small = shtools_to_healpy(clm, lmax=5)
    big = shtools_to_healpy(clm, lmax=LMAX + 4)
    l, m = hp.Alm.getlm(5)
    np.testing.assert_allclose(small, full[hp.Alm.getidx(LMAX, l, m)])
    l, m = hp.Alm.getlm(LMAX + 4)
    np.testing.assert_allclose(big[l <= LMAX], full[hp.Alm.getidx(LMAX, l[l <= LMAX], m[l <= LMAX])])
    assert np.all(big[l > LMAX] == 0)


# ----------------------------------------------------------------------------------
# 4. End-to-end: grid -> pyshtools expand -> healpy, vs healpy map2alm
# ----------------------------------------------------------------------------------
def test_end_to_end_against_map2alm():
    """End-to-end: DH grid -> pyshtools expand -> healpy agrees with healpy map2alm.

    Mirrors the AEP workflow: sample a band-limited field on a DH grid, expand with
    pyshtools (real, ortho, csphase=-1), convert, and compare with healpy map2alm
    of the same field sampled on a HEALPix grid.
    """
    lmax, nside = 16, 64
    truth = random_real_clm(lmax, '4pi', 1, seed=7)
    grid = truth.expand(grid='DH2', extend=True)
    clm = grid.expand(normalization='ortho', csphase=-1, lmax_calc=lmax)
    alm_conv = shtools_to_healpy(clm)

    theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    hmap = truth.expand(lat=90 - np.degrees(theta), lon=np.degrees(phi), degrees=True)
    alm_hp = hp.map2alm(hmap, lmax=lmax, iter=10)
    np.testing.assert_allclose(alm_conv, alm_hp, atol=1e-6 * np.abs(alm_hp).max())


# ----------------------------------------------------------------------------------
# 5. The pitfalls in the original script, documented as tests
# ----------------------------------------------------------------------------------
def test_pitfall_naive_real_to_complex_is_wrong():
    """Using alm = C + iS is wrong: it needs (C - iS)/sqrt(2) for m > 0."""
    clm = random_real_clm(LMAX, 'ortho', -1)
    l, m = hp.Alm.getlm(LMAX)
    naive = clm.coeffs[0][l, m] + 1j * clm.coeffs[1][l, m]
    assert not np.allclose(naive, shtools_to_healpy(clm))
    fixed = np.where(m == 0, naive, np.conj(naive) / np.sqrt(2))
    np.testing.assert_allclose(fixed, shtools_to_healpy(clm), atol=1e-14)


def test_pitfall_makegridpoint_default_csphase(pix):
    """MakeGridPoint defaults to csphase=1: omitting csphase=-1 rotates the field by 180 deg.

    Evaluating csphase=-1 coefficients without passing csphase=-1 flips the sign of every
    odd-m term, i.e. rotates the field by 180 deg in azimuth.
    """
    _, phi, lat, lon = pix
    clm = random_real_clm(LMAX, 'ortho', -1)
    right = pysh.expand.MakeGridPoint(clm.coeffs, lat, lon, norm=4, csphase=-1)
    wrong = pysh.expand.MakeGridPoint(clm.coeffs, lat, lon, norm=4)
    assert not np.allclose(right, wrong)
    rotated = pysh.expand.MakeGridPoint(clm.coeffs, lat, (lon + 180) % 360, norm=4, csphase=-1)
    np.testing.assert_allclose(wrong, rotated, atol=1e-10)
    np.testing.assert_allclose(right, hp.alm2map(shtools_to_healpy(clm), NSIDE, lmax=LMAX), atol=1e-10)


# ----------------------------------------------------------------------------------
# 6. healpy-free alm index helpers, and .sph data
# ----------------------------------------------------------------------------------
@pytest.mark.parametrize('lmax', [0, 1, 7, 50])
def test_alm_index_helpers_match_healpy(lmax):
    """The index helpers reproduce healpy's Alm helpers (mmax = lmax)."""
    l, m = alm_getlm(lmax)
    lh, mh = hp.Alm.getlm(lmax)
    np.testing.assert_array_equal(l, lh)
    np.testing.assert_array_equal(m, mh)
    assert alm_getsize(lmax) == hp.Alm.getsize(lmax)
    assert alm_getlmax(alm_getsize(lmax)) == lmax
    np.testing.assert_array_equal(alm_getidx(lmax, l, m), hp.Alm.getidx(lmax, l, m))
    assert alm_getlmax(alm_getsize(lmax) + 1) == -1


def test_sph_potentials_to_healpy(pix):
    """The .sph potentials U, V convert to healpy alm pairs and back without changing the far field.

    U and V are complex fields, so each needs two healpy alm arrays (real and imaginary parts).
    """
    _, _, lat, lon = pix
    d = read_sph_shcoeffs('tests/sph/hertzian_z_dip_array_FarField1_299MHz.sph')
    lmax = d['nmax']
    out = {}
    for key in ('U', 'V'):
        clm = d[key]
        re, im = shtools_to_healpy(clm, complex_field=True)
        f_hp = hp.alm2map(re, NSIDE, lmax=lmax) + 1j * hp.alm2map(im, NSIDE, lmax=lmax)
        f_sh = clm.expand(lat=lat, lon=lon, degrees=True)
        np.testing.assert_allclose(f_hp, f_sh, atol=1e-10 * np.abs(f_sh).max())
        out[key] = healpy_to_shtools(re, alm_imag=im, normalization='ortho', csphase=-1, kind='complex')
    theta, phi = np.deg2rad(np.arange(0, 181, 5)), np.deg2rad(np.arange(0, 360, 5))
    E0 = shcoeffs_to_farfield(d['U'], d['V'], theta, phi)
    E1 = shcoeffs_to_farfield(out['U'], out['V'], theta, phi)
    for a, b in zip(E1, E0, strict=True):
        np.testing.assert_allclose(a, b, atol=1e-10 * np.abs(b).max())
