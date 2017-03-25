def alm_normal(lmax, mean=0.0, sigma=1.0):
    from numpy.random import normal

    a_coef = normal(mean, sigma, size=(lmax + 1, lmax + 1))
    b_coef = normal(mean, sigma, size=(lmax + 1, lmax + 1))

    for m in xrange(0, lmax + 1):
        for l in xrange(0, m):
            a_coef[m][l] = 0.0
            b_coef[m][l] = 0.0
    for l in xrange(0, lmax + 1):
        b_coef[0][l] = 0.0

    a_coef[0][0] = 0.0
    b_coef[0][0] = 0.0
    a_coef[0][1] = 0.0
    a_coef[1][1] = 0.0
    b_coef[0][1] = 0.0
    b_coef[1][1] = 0.0

    return a_coef, b_coef


def alm_l(lmax, per_sigma=0.2, k0=1.0, k1=1.0):
    from numpy.random import normal
    from numpy import zeros

    a_coef = zeros((lmax + 1, lmax + 1))
    b_coef = zeros((lmax + 1, lmax + 1))

    for l in xrange(0, lmax + 1):

        mean = 1/(1/k0 + l)**k1
        sigma = per_sigma * mean

        for m in xrange(0, l + 1):
            a_coef[m][l] = normal(mean, sigma)
            b_coef[m][l] = normal(mean, sigma)

    for m in xrange(0, lmax + 1):
        for l in xrange(0, m):
            a_coef[m][l] = 0.0
            b_coef[m][l] = 0.0
    for l in xrange(0, lmax + 1):
        b_coef[0][l] = 0.0

    a_coef[0][0] = 0.0
    b_coef[0][0] = 0.0
    a_coef[0][1] = 0.0
    a_coef[1][1] = 0.0
    b_coef[0][1] = 0.0
    b_coef[1][1] = 0.0

    return a_coef, b_coef


def cl(a, b, lmax):
    from numpy import zeros

    c = zeros((lmax + 1))

    for l in xrange(0, lmax + 1):

        c_sum = 0.0

        for m in xrange(0, l + 1):
            c_sum = c_sum + a[m][l] * a[m][l] + b[m][l] * b[m][l]

        c[l] = c_sum / (2.0 * l + 1.0)

    return c


def correlations_from_map(x, y, field):
    from math import cos, sqrt

    n = field.shape[0] - 1

    a = 0.0
    na = 0.0

    for i in xrange(0, n):
        for j in xrange(1, n / 2):
            a += cos(y[i][j]) * field[i][j] * field[i][j]
            na += cos(y[i][j])

    return sqrt(a / na)


def correlations_from_alm(a, b, lmax):
    from math import sqrt, pi

    c = cl(a, b, lmax)

    sigma_0 = 0.0
    for l in xrange(0, lmax + 1):
        sigma_0 += (2.0 * l + 1.0) * c[l]
    sigma_0 = sqrt(sigma_0 / 4.0 / pi)

    sigma_1 = 0.0
    for l in xrange(0, lmax + 1):
        sigma_1 += l * (l + 1.0) * (2.0 * l + 1.0) * c[l]
    sigma_1 = sqrt(sigma_1 / 4.0 / pi)

    sigma_2 = 0.0
    for l in xrange(0, lmax + 1):
        sigma_2 += (l + 2.0) * (l - 1.0) * l * (l + 1.0) * (2.0 * l + 1.0) * c[l]
    sigma_2 = sqrt(sigma_2 / 4.0 / pi)

    return sigma_0, sigma_1, sigma_2


def plot_al_m(alm_normal_format, l_max):
    from numpy import linspace, real
    from pylab import plot, show

    i = linspace(0, l_max, l_max)
    plot(i, real(alm_normal_format[1]))

    show()


def print_al_m_c(alm_normal_format, l_max):
    from numpy import linspace, real
    from pylab import plot, show

    i = linspace(0, l_max, l_max)
    plot(i, real(alm_normal_format[1]) * real(alm_normal_format[1]) * i * (i + 1))

    show()


def plot_bl_m(alm_normal_format, l_max):
    from numpy import linspace, imag
    from pylab import plot, show

    i = linspace(0, l_max, l_max)
    plot(i, imag(alm_normal_format[1]))

    show()


def print_bl_m_c(alm_normal_format, l_max):
    from numpy import linspace, imag
    from pylab import plot, show

    i = linspace(0, l_max, l_max)
    plot(i, imag(alm_normal_format[1]) * imag(alm_normal_format[1]) * i * (i + 1))

    show()
