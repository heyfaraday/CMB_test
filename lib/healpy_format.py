def healpy_alm(filename, lmax):
    import healpy as hp
    import numpy as np

    alms_from_file = hp.fitsfunc.read_alm(filename)
    lvals, mvals = hp.sphtfunc.Alm.getlm(lmax, np.arange(hp.sphtfunc.Alm.getsize(lmax)))

    alm_normal = np.zeros((lmax + 1, lmax + 1), dtype=complex)  # It's out formating

    for i in xrange(0, np.size(lvals)):

        m = mvals[i]
        l = lvals[i]
        alm_normal[m][l] = alms_from_file[i]

    return alm_normal
