import matplotlib.pyplot as plt

import numpy as np
import healpy as hp

map_I = hp.read_map('data/COM_CMB_IQU-smica_1024_R2.02_full.fits')
hp.mollview(map_I, norm='hist', min=-0.1, max=0.1, xsize=2000)
plt.show()

map_Q = hp.read_map('data/COM_CMB_IQU-smica_1024_R2.02_full.fits', field=1)
hp.mollview(map_Q, norm='hist', min=-0.01, max=0.01, xsize=2000)
plt.show()

map_U = hp.read_map('data/COM_CMB_IQU-smica_1024_R2.02_full.fits', field=2)
hp.mollview(map_U, norm='hist', min=-0.01, max=0.01, xsize=2000)
plt.show()

cl_I = hp.anafast(map_I, lmax=2048)
plt.show()

cl_Q = hp.anafast(map_Q, lmax=2048)
plt.show()

cl_U = hp.anafast(map_U, lmax=2048)
plt.show()

ell = np.arange(len(cl_I))

plt.figure(figsize=(5, 5))
plt.plot(ell, ell * (ell + 1) * cl_I)
plt.xlabel('ell')
plt.ylabel('ell(ell+1)cl_I')
plt.grid()
plt.show()

plt.figure(figsize=(5, 5))
plt.plot(ell, ell * (ell + 1) * cl_Q)
plt.xlabel('ell')
plt.ylabel('ell(ell+1)cl_Q')
plt.grid()
plt.show()

plt.figure(figsize=(5, 5))
plt.plot(ell, ell * (ell + 1) * cl_U)
plt.xlabel('ell')
plt.ylabel('ell(ell+1)cl_U')
plt.grid()
plt.show()
