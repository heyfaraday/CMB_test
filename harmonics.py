import healpy as hp
import numpy as np

lmax = 512

map_I = hp.read_map('data/COM_CMB_IQU-smica_1024_R2.02_full.fits', field=1)

alm_I = hp.sphtfunc.map2alm(map_I, lmax=lmax)

hp.fitsfunc.write_alm('planck_2_dir/planck_2_Q_512.fits', alm_I)
file_map = open('planck_2_dir/planck_2_Q_norm_512.dat', 'w')

alms_from_file = hp.fitsfunc.read_alm('planck_2_dir/planck_2_Q_512.fits')
lvals, mvals = hp.sphtfunc.Alm.getlm(lmax, np.arange(hp.sphtfunc.Alm.getsize(lmax)))

alm_normal = np.zeros((lmax + 1, lmax + 1), dtype=complex)

for i in xrange(0, np.size(lvals)):
    m = mvals[i]
    l = lvals[i]
    alm_normal[m][l] = alms_from_file[i]

for m in xrange(0, lmax+1):
    for l in xrange(0, lmax+1):
        file_map.write(repr(alm_normal[m][l]) + '    ' + repr(m) + '    ' + repr(l) + '\n')
