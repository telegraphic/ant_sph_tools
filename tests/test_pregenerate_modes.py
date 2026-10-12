import numpy as np
import pylab as plt
from scipy.io import loadmat
import time

from ant_sph_tools import read_sph
from ant_sph_tools.pregenerate_modes import pregenerate

sph_data = read_sph('tests/sph-aavs/AAVS1d5_single_element_rev1_110MHz_FarField1(0).sph')
theta = np.deg2rad(np.arange(0, 181)/2)
phi = np.deg2rad(np.arange(0, 721)/2)

t0 = time.time()
E_th, E_ph, mode_counter = pregenerate(sph_data['Q'], 25, 25, theta, phi)
t1 = time.time()

print(f"{t1-t0:.2f}s")
