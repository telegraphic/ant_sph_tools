# Spherical harmonics in pyshtools and healpy

[SHTOOLS](https://shtools.github.io/SHTOOLS/) (`pyshtools`) and [healpy](https://healpy.readthedocs.io/) both provide spherical harmonic tools, but they use different normalizations, phase conventions, real/complex bases and storage layouts, so coefficients cannot be copied from one to the other directly. This page summarises the differences and documents `ant_sph_tools.healpix`, which converts between them. Every statement here is checked numerically in `tests/test_healpix.py`.

## Summary

| | pyshtools | healpy |
|---|---|---|
| Harmonics | real (default) or complex, $Y_{\ell m}$ | complex, $Y_{\ell m}$ |
| Coefficients | real $C_{\ell m}, S_{\ell m}$; complex $c_{\ell m}$ | $a_{\ell m}$ |
| Normalization | `'4pi'` (default), `'ortho'`, `'schmidt'`, `'unnorm'` | orthonormal |
| Condon-Shortley phase $(-1)^m$ | excluded by default (`csphase=1`) | **included** |
| Field | real or complex | real |
| Storage | `(2, lmax+1, lmax+1)` array | 1-D, $m \ge 0$ only, ordered by $m$ then $\ell$ |
| Positions | latitude, longitude (degrees) | colatitude $\theta$, longitude $\phi$ (radians) |

In both, $\theta$ is colatitude (0 at the North pole) and $\phi$ is longitude, increasing eastward. Below, $P_{\ell m}$ is the associated Legendre function *without* the Condon-Shortley phase; where the phase is used it is written explicitly, and $[\mathrm{cs}]$ marks a factor that is present only for `csphase=-1`.

> **Key fact:** pyshtools complex harmonics with `normalization='ortho'` and `csphase=-1` are identical to the healpy $Y_{\ell m}$ (for both positive and negative $m$). In that convention $a_{\ell m} = c_{\ell m}$, i.e. `alm = clm.coeffs[0][l, m]`.

## Definitions

### pyshtools, real harmonics (`kind='real'`, the default)

```math
f(\theta,\phi) = \sum_{\ell=0}^{\infty}\sum_{m=0}^{\ell}\big[C_{\ell m}\cos m\phi + S_{\ell m}\sin m\phi\big]\,\bar P_{\ell m}(\cos\theta),
\qquad
\bar P_{\ell m}(\mu) = \sqrt{(2-\delta_{m0})(2\ell+1)\frac{(\ell-m)!}{(\ell+m)!}}\;(-1)^{m\,[\mathrm{cs}]}\,P_{\ell m}(\mu)
```

* $4\pi$ normalized: $\int Y_{\ell m}^2\,d\Omega = 4\pi$.
* Note the $\sqrt{2-\delta_{m0}}$: real harmonics with $m>0$ are $\sqrt 2$ larger than the complex ones.
* Storage: `coeffs[0, l, m]` $= C_{\ell m}$, `coeffs[1, l, m]` $= S_{\ell m}$.

### pyshtools, complex harmonics (`kind='complex'`)

```math
f(\theta,\phi) = \sum_{\ell=0}^{\infty}\sum_{m=-\ell}^{\ell} c_{\ell m} Y_{\ell m}(\theta,\phi),
\qquad
Y_{\ell m} = \bar P_{\ell m}(\cos\theta)\,e^{im\phi},
\qquad
\bar P_{\ell m}(\mu) = \sqrt{(2\ell+1)\frac{(\ell-m)!}{(\ell+m)!}}\;(-1)^{m\,[\mathrm{cs}]}\,P_{\ell m}(\mu)
```

```math
\bar P_{\ell,-m} = (-1)^m \bar P_{\ell m}, \qquad Y_{\ell,-m} = (-1)^m Y^*_{\ell m}
```

* $4\pi$ normalized, *without* the $\sqrt{2-\delta_{m0}}$.
* Storage: `coeffs[0, l, m]` $= c_{\ell m}$, `coeffs[1, l, m]` $= c_{\ell,-m}$.
* A real field has $c_{\ell,-m} = (-1)^m c^*_{\ell m}$.

### healpy

```math
f(\theta,\phi) = \sum_{\ell=0}^{\ell_{\max}}\sum_{m=-\ell}^{\ell} a_{\ell m} Y_{\ell m}(\theta,\phi),
\qquad
Y_{\ell m} = \lambda_{\ell m}(\cos\theta)\,e^{im\phi},
\qquad
\lambda_{\ell m}(x) = (-1)^m\sqrt{\frac{2\ell+1}{4\pi}\frac{(\ell-m)!}{(\ell+m)!}}\,P_{\ell m}(x)
```

```math
Y_{\ell,-m} = (-1)^m Y^*_{\ell m}, \qquad a_{\ell,-m} = (-1)^m a^*_{\ell m}
```

* Orthonormal: $\int |Y_{\ell m}|^2\,d\Omega = 1$.
* The Condon-Shortley phase is **included** (it does *not* match the pyshtools default). For example, $a_{11}=1$ gives the map $2\,\mathrm{Re}(Y_{11}) = -\sqrt{3/2\pi}\,\sin\theta\cos\phi$.
* Maps are real, so only $m\ge 0$ is stored, at index $m(2\ell_{\max}+1-m)/2 + \ell$ (`hp.Alm.getidx`, or `ant_sph_tools.healpix.alm_getidx`).

## Converting pyshtools to healpy

Four independent differences each break a naive `alm = coeffs[0][l, m]`:

1. **Phase** — pyshtools excludes the Condon-Shortley phase by default; healpy includes it. Odd-$m$ coefficients flip sign.
2. **Normalization** — pyshtools defaults to $4\pi$; healpy is orthonormal. Coefficients scale by $F_{\ell m}$ (table below).
3. **Real vs complex** — pyshtools defaults to *real* harmonics, which split each order into cosine and sine terms *and* carry an extra $\sqrt 2$ for $m>0$. The complex coefficient is $(C_{\ell m} - iS_{\ell m})/\sqrt 2$, **not** $C_{\ell m} + iS_{\ell m}$. The wrong sign mirrors the map in longitude ($\phi\to-\phi$); the missing $1/\sqrt 2$ over-weights all $m>0$ terms relative to $m=0$ (which can, e.g., drive a positive-definite beam pattern negative).
4. **Real fields only** — one healpy alm array describes a real map. A complex field $f = g + ih$ needs two arrays, one each for $g$ and $h$.

For a real field:

```math
a_{\ell m} = s_m F_{\ell m} \times
\begin{cases}
c_{\ell m} & \text{complex pyshtools coefficients}\\
C_{\ell 0} & \text{real pyshtools coefficients, } m=0\\
(C_{\ell m} - iS_{\ell m})/\sqrt 2 & \text{real pyshtools coefficients, } m>0
\end{cases}
```

where $s_m = (-1)^m$ for pyshtools `csphase=1` (the default) and $s_m = 1$ for `csphase=-1`. The inverse, for real pyshtools output, is $C_{\ell m} = \sqrt 2\,\mathrm{Re}(a_{\ell m})/(s_m F_{\ell m})$, $S_{\ell m} = -\sqrt 2\,\mathrm{Im}(a_{\ell m})/(s_m F_{\ell m})$ for $m>0$, and $C_{\ell 0} = a_{\ell 0}/F_{\ell 0}$.

| pyshtools normalization | $F_{\ell m}$ (real) | $F_{\ell m}$ (complex) |
|---|---|---|
| `'ortho'` | $1$ | $1$ |
| `'4pi'` (default) | $\sqrt{4\pi}$ | $\sqrt{4\pi}$ |
| `'schmidt'` | $\sqrt{4\pi/(2\ell+1)}$ | $\sqrt{4\pi/(2\ell+1)}$ |
| `'unnorm'` | $\sqrt{\frac{4\pi}{(2-\delta_{m0})(2\ell+1)}\frac{(\ell+m)!}{(\ell-m)!}}$ | $\sqrt{\frac{4\pi}{2\ell+1}\frac{(\ell+m)!}{(\ell-m)!}}$ |

$F_{\ell m}$ converts pyshtools coefficients to orthonormal coefficients of the same kind (`norm_to_ortho_factor`). For real harmonics the $\sqrt{2-\delta_{m0}}$ is already inside $\bar P_{\ell m}$ for all but the unnormalized case, and is accounted for by the $1/\sqrt 2$ above.

For a complex field, with $a^{\pm}_{\ell m} = s_m F_{\ell m}\,c_{\ell,\pm m}$ (and $a^-_{\ell 0} = a^+_{\ell 0}$), the healpy coefficients of the real and imaginary parts are

```math
a^{(g)}_{\ell m} = \tfrac12\big[a^+_{\ell m} + (-1)^m a^{-*}_{\ell m}\big],
\qquad
a^{(h)}_{\ell m} = \tfrac{1}{2i}\big[a^+_{\ell m} - (-1)^m a^{-*}_{\ell m}\big],
```

so that `f = hp.alm2map(alm_g, nside) + 1j * hp.alm2map(alm_h, nside)`.

## Usage

```python
import healpy as hp
import pyshtools as pysh
from ant_sph_tools import healpy_to_shtools, shtools_to_healpy

clm = pysh.SHGrid.from_array(data).expand()          # pyshtools defaults: 4pi, csphase=1, real
alm = shtools_to_healpy(clm)                          # healpy convention
hmap = hp.alm2map(alm, nside, lmax=clm.lmax)          # same field on a HEALPix grid

alm = shtools_to_healpy(coeffs, normalization='ortho', csphase=-1)   # raw (2, L+1, L+1) array
alm_re, alm_im = shtools_to_healpy(clm_complex, complex_field=True)   # complex-valued field
clm2 = healpy_to_shtools(alm, normalization='4pi', csphase=1, kind='real')
```

Without `ant_sph_tools`, the shortest correct route is to let pyshtools change convention and then reindex:

```python
l, m = hp.Alm.getlm(clm.lmax)
alm = clm.convert(normalization='ortho', csphase=-1, kind='complex').coeffs[0][l, m]
```

`ant_sph_tools.healpix` does not import healpy; `alm_getlm`, `alm_getidx`, `alm_getsize` and `alm_getlmax` reproduce the `hp.Alm` helpers (for `mmax == lmax`).

### Spherical mode (.sph) data

`read_sph_shcoeffs()` returns the far field of a .sph file as two pyshtools potentials, `U` (TM modes) and `V` (TE modes), with $\mathbf E = \nabla_s U - \hat r \times \nabla_s V$ (see the README). By default these use `'ortho'`, `csphase=-1`, i.e. the healpy basis. They are complex fields, so each converts to a pair of healpy alm arrays:

```python
from ant_sph_tools import read_sph_shcoeffs, shtools_to_healpy

d = read_sph_shcoeffs('tests/sph/dipole_FarField1_299MHz.sph')
U_re, U_im = shtools_to_healpy(d['U'], complex_field=True)
V_re, V_im = shtools_to_healpy(d['V'], complex_field=True)
```

## Common pitfalls

* `pysh.expand.MakeGridPoint` (and `MakeGridPointC`) default to `csphase=1`. Evaluating `csphase=-1` coefficients without passing `csphase=-1` flips every odd-$m$ term, which rotates the field by 180° in azimuth. For patterns that are nearly symmetric under that rotation (e.g. dipole power patterns) this is easy to miss.
* `alm = C + 1j * S` from real pyshtools coefficients is wrong in both sign and scale; use $(C - iS)/\sqrt 2$ for $m>0$ (item 3 above).
* Comparing against `hp.map2alm` is only accurate to its quadrature error; use `iter` > 0 (and a sufficiently high `nside`) when using it as a reference.
