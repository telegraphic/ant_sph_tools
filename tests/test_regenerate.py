import time

import numpy as np
import pylab as plt
from scipy.io import loadmat

from ant_sph_tools import read_sph
from ant_sph_tools.farfield_pregen import regenerate

sph_data = read_sph('tests/sph-aavs/AAVS1d5_single_element_rev1_110MHz_FarField1(0).sph')
theta = np.deg2rad(np.arange(0, 181)/2)
phi = np.deg2rad(np.arange(0, 721)/2)

t0 = time.time()
E_th, E_ph, mode_counter = regenerate(sph_data['Q'], 20, 20, theta, phi)
t1 = time.time()

print(f"{t1-t0:.2f}s")

mat = loadmat('tests/sph-aavs/AAVS1d5_single_element_110MHz_Xpol.mat')

plt.figure(figsize=(8,8))

plt.subplot(2,2,1)
plt.imshow(np.abs(E_ph), aspect='auto')
plt.title("Python $E_{\\phi}$")
plt.colorbar()

plt.subplot(2,2,2)
plt.imshow(np.abs(mat['Ephi']), aspect='auto')
plt.title("MATLAB $E_{\\phi}$")
plt.colorbar()

plt.subplot(2,2,3)
plt.imshow(np.abs(E_th), aspect='auto')
plt.title("Python $E_{\\theta}$")
plt.colorbar()

plt.subplot(2,2,4)
plt.imshow(np.abs(mat['Etheta']), aspect='auto')
plt.title("MATLAB $E_{\\theta}$")
plt.colorbar()
plt.suptitle("AAVS1d5 single_element 110MHz - X pol")
plt.tight_layout()
plt.show()

assert np.allclose(E_ph, mat['Ephi']), np.allclose(E_th, mat['Etheta'])