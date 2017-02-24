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


def inverse_healpix(size, f, p_, back_f, back_fa, back_fb, l_max_inv):
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
