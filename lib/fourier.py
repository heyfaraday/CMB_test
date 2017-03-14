def direct_f(size, p_, fa, fb, cos_coef, sin_coef, l_max_dir, diff=False):
    from numpy import zeros, real, imag
    from numpy.fft import fft

    f = zeros((size + 1, size / 2 + 1))

    func1 = 0.0
    func2 = 0.0

    if diff:

        for j in xrange(1, size / 2):

            for m in xrange(0, l_max_dir + 1):
                for l in xrange(m, l_max_dir + 1):
                    func1 += cos_coef[m][l] * p_[j][m][l]
                    func2 += sin_coef[m][l] * p_[j][m][l]

                fa[j][m] = func1
                fb[j][m] = func2

                func1 = 0.0
                func2 = 0.0

            t = - imag(fft(fa[j])) + real(fft(fb[j]))

            f[0:size, j] = t[:]
            f[size][j] = f[0][j]

    else:

        for j in xrange(1, size / 2):

            for m in xrange(0, l_max_dir + 1):
                for l in xrange(m, l_max_dir + 1):
                    func1 += cos_coef[m][l] * p_[j][m][l]
                    func2 += sin_coef[m][l] * p_[j][m][l]

                fa[j][m] = func1
                fb[j][m] = func2

                func1 = 0.0
                func2 = 0.0

            t = real(fft(fa[j])) + imag(fft(fb[j]))

            f[0:size, j] = t[:]
            f[size][j] = f[0][j]

    return f


def inverse_f(size, f, p_, back_f, back_fa, back_fb, l_max_inv):
    from numpy import zeros, real, imag
    from numpy.fft import ifft
    from math import pi, sin

    f = f.T

    back_cos_coef = zeros((l_max_inv + 1, l_max_inv + 1))
    back_sin_coef = zeros((l_max_inv + 1, l_max_inv + 1))

    norm = 0.0
    for j in xrange(1, size / 2):
        theta = 2.0 * pi * j / size
        norm += sin(theta)

    for my_m in xrange(0, l_max_inv + 1):
        for my_l in xrange(my_m, l_max_inv + 1):

            sum_a = 0.0
            sum_b = 0.0

            for j in xrange(1, size / 2):
                theta = 2.0 * pi * j / size

                back_f[j] = ifft(f[j][0:size])

                back_fa[j] = real(back_f[j])
                back_fb[j] = imag(back_f[j])

                sum_a += p_[j][my_m][my_l] * back_fa[j][my_m] * sin(theta) / norm * 4 * pi
                sum_b -= p_[j][my_m][my_l] * back_fb[j][my_m] * sin(theta) / norm * 4 * pi

            back_cos_coef[my_m][my_l] = sum_a
            back_sin_coef[my_m][my_l] = sum_b

    return back_cos_coef, back_sin_coef


def coef_1(in_l, in_m):
    from math import sqrt

    if in_l != 0:
        return sqrt((in_l - in_m) * (2.0 * in_l + 1.0)
                    / ((in_l + in_m) * (2.0 * in_l - 1.0)))
    if in_l == 0:
        return 0.0


def legendre(j, n, l_max):
    from math import pi, sqrt, sin, cos
    from numpy import zeros

    theta = 2.0 * pi * j / n

    p_ = zeros((l_max + 1, l_max + 1))

    p_[0][0] = 1.0 / sqrt(4.0 * pi)

    for m in xrange(0, l_max):
        p_[m + 1][m + 1] = - p_[m][m] * sin(theta) * sqrt(2.0 * m + 3.0) / sqrt(2.0 * m + 2.0)

    for m in xrange(0, l_max):
        p_[m][m + 1] = p_[m][m] * cos(theta) * sqrt(2.0 * m + 3.0)

    for m in xrange(0, l_max - 1):
        for l in xrange(m + 2, l_max + 1):
            p_[m][l] = ((2.0 * l - 1.0) * sqrt((l - m) * (2.0 * l + 1.0)) / sqrt((l + m) * (2.0 * l - 1.0)) *
                           p_[m][l - 1] * cos(theta) - (l + m - 1.0) * sqrt((l - m) * (l - 1.0 - m) *
                                                                               (2.0 * l + 1.0)) / sqrt(
                (l + m) * (l - 1.0 + m) * (2.0 * l - 3.0)) * p_[m][l - 2]) / \
                          (l - m)

    for m in xrange(1, l_max + 1):
        for l in xrange(m, l_max + 1):
            p_[m][l] *= sqrt(2.0)

    return p_


def legendre_x(j, n, l_max):
    from math import pi, sin
    from numpy import zeros

    theta = 2.0 * pi * j / n

    f_x = zeros((l_max + 1, l_max + 1))
    p_ = legendre(j, n, l_max)

    for l in xrange(2, l_max + 1):
        for m in xrange(0, l + 1):
            f_x[m][l] = m * p_[m][l] / sin(theta)
    return f_x


def legendre_y(j, n, l_max):
    from math import pi, sin, cos
    from numpy import zeros

    theta = 2.0 * pi * j / n

    f_y = zeros((l_max + 1, l_max + 1))
    p_ = legendre(j, n, l_max)

    for l in xrange(2, l_max + 1):
        for m in xrange(0, l + 1):
            f_y[m][l] = l * cos(theta) / sin(theta) * p_[m][l] - \
                        1.0 / sin(theta) * (l + m) * coef_1(l, m) * p_[m][l - 1]
    return f_y


def direct_f_int(size, fa, fb, cos_coef, sin_coef, l_max_dir, sign, diff=False):
    from numpy import zeros, real, imag
    from numpy.fft import fft

    f = zeros((size + 1, size / 2 + 1))

    func1 = 0.0
    func2 = 0.0

    if diff:

        for j in xrange(1, size / 2):

            if sign[0] == '0':
                p_ = legendre(j, size, l_max_dir)
            elif sign[0] == 'x':
                p_ = legendre_x(j, size, l_max_dir)
            elif sign[0] == 'y':
                p_ = legendre_y(j, size, l_max_dir)
            else:
                print 'Not valid'

            for m in xrange(0, l_max_dir + 1):
                for l in xrange(m, l_max_dir + 1):
                    func1 += cos_coef[m][l] * p_[m][l]
                    func2 += sin_coef[m][l] * p_[m][l]

                fa[j][m] = func1
                fb[j][m] = func2

                func1 = 0.0
                func2 = 0.0

            t = - imag(fft(fa[j])) + real(fft(fb[j]))

            f[0:size, j] = t[:]
            f[size][j] = f[0][j]

    else:

        for j in xrange(1, size / 2):

            if sign[0] == '0':
                p_ = legendre(j, size, l_max_dir)
            elif sign[0] == 'x':
                p_ = legendre_x(j, size, l_max_dir)
            elif sign[0] == 'y':
                p_ = legendre_y(j, size, l_max_dir)
            else:
                print 'Not valid'

            for m in xrange(0, l_max_dir + 1):
                for l in xrange(m, l_max_dir + 1):
                    func1 += cos_coef[m][l] * p_[m][l]
                    func2 += sin_coef[m][l] * p_[m][l]

                fa[j][m] = func1
                fb[j][m] = func2

                func1 = 0.0
                func2 = 0.0

            t = real(fft(fa[j])) + imag(fft(fb[j]))

            f[0:size, j] = t[:]
            f[size][j] = f[0][j]

    return f
