from math import sqrt, pi, sin, cos

import numpy as np

from lib.fourier import direct_f, inverse_f

from lib import cmbplot

L_max_field = 4
L_max_polynom = 4
L_max_back = 4
N = 32


def coef_1(in_l, in_m):
    if in_l != 0:
        return sqrt((in_l - in_m) * (2.0 * in_l + 1.0)
                    / ((in_l + in_m) * (2.0 * in_l - 1.0)))
    if in_l == 0:
        return 0.0


# P_ generation
P_ = np.zeros((N / 2 + 1, L_max_polynom + 1, L_max_polynom + 1))
for j in xrange(1, N / 2):

    theta = 2.0 * pi * j / N

    P_[j][0][0] = 1.0 / sqrt(4.0 * pi)

    for m in xrange(0, L_max_polynom):
        P_[j][m + 1][m + 1] = - P_[j][m][m] * sin(theta) * sqrt(2.0 * m + 3.0) / sqrt(2.0 * m + 2.0)

    for m in xrange(0, L_max_polynom):
        P_[j][m][m + 1] = P_[j][m][m] * cos(theta) * sqrt(2.0 * m + 3.0)

    for m in xrange(0, L_max_polynom - 1):
        for l in xrange(m + 2, L_max_polynom + 1):
            P_[j][m][l] = ((2.0 * l - 1.0) * sqrt((l - m) * (2.0 * l + 1.0)) / sqrt((l + m) * (2.0 * l - 1.0)) *
                            P_[j][m][l - 1] * cos(theta) - (l + m - 1.0) * sqrt((l - m) * (l - 1.0 - m) *
                            (2.0 * l + 1.0)) / sqrt((l + m) * (l - 1.0 + m) * (2.0 * l - 3.0)) * P_[j][m][l - 2]) / \
                            (l - m)

    for m in xrange(1, L_max_polynom + 1):
        for l in xrange(m, L_max_polynom + 1):
            P_[j][m][l] *= sqrt(2.0)

# F_x generation - np.imag + np.real
F_x = np.zeros((N / 2 + 1, L_max_polynom + 1, L_max_polynom + 1))
for j in xrange(1, N / 2):

    theta = 2.0 * pi * j / N

    for l in xrange(2, L_max_polynom + 1):
        for m in xrange(0, l + 1):
            F_x[j][m][l] = m * P_[j][m][l] / sin(theta)

# F_y generation - np.real + np.imag
F_y = np.zeros((N / 2 + 1, L_max_polynom + 1, L_max_polynom + 1))
for j in xrange(1, N / 2):

    theta = 2.0 * pi * j / N

    for l in xrange(2, L_max_polynom + 1):
        for m in xrange(0, l + 1):
            F_y[j][m][l] = l * cos(theta) / sin(theta) * P_[j][m][l] - \
                           1.0 / sin(theta) * (l + m) * coef_1(l, m) * P_[j][m][l - 1]

# F_xy generation - np.imag + np.real
F_xy = np.zeros((N / 2 + 1, L_max_polynom + 1, L_max_polynom + 1))
for j in xrange(1, N / 2):

    theta = 2.0 * pi * j / N

    for l in xrange(2, L_max_polynom + 1):
        for m in xrange(0, l + 1):
            F_xy[j][m][l] = m / sin(theta) * ((1.0 / sin(theta)) * (l + m) * P_[j][m][l - 1] * coef_1(l, m) -
                                                       (l - 1.0) * cos(theta) / sin(theta) * P_[j][m][l])

# F_xx_1 generation - np.real + np.real
F_xx_1 = np.zeros((N / 2 + 1, L_max_polynom + 1, L_max_polynom + 1))
for j in xrange(1, N / 2):

    theta = 2.0 * pi * j / N

    for l in xrange(2, L_max_polynom + 1):
        for m in xrange(0, l + 1):
            F_xx_1[j][m][l] = - m * m * P_[j][m][l] / (sin(theta) * sin(theta))

# F_xx_2 generation - np.real + np.real
F_xx_2 = np.zeros((N / 2 + 1, L_max_polynom + 1, L_max_polynom + 1))
for j in xrange(1, N / 2):

    theta = 2.0 * pi * j / N

    for l in xrange(2, L_max_polynom + 1):
        for m in xrange(0, l + 1):
            F_xx_2[j][m][l] = (l * cos(theta) / sin(theta) * P_[j][m][l] - 1.0 / sin(theta) * (l + m) * coef_1(l, m) *
                               P_[j][m][l - 1]) * cos(theta) / sin(theta)

# F_yy generation - np.real + np.real
F_yy = np.zeros((N / 2 + 1, L_max_polynom + 1, L_max_polynom + 1))
for j in xrange(1, N / 2):

    theta = 2.0 * pi * j / N

    for l in xrange(2, L_max_polynom + 1):
        for m in xrange(0, l + 1):
            F_yy[j][m][l] = 0.5 / sin(theta) * ((1.0 / sin(theta)) * (l * l * cos(2.0 * theta) -
                            (l + 2.0) * l + 2.0 * m * m) * P_[j][m][l] + 2.0 * (l + m) * cos(theta) /
                            sin(theta) * coef_1(l, m) * P_[j][m][l - 1])

x = np.zeros((N + 1, N / 2 + 1))
y = np.zeros((N + 1, N / 2 + 1))

for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        x[i][j] = (2.0 * i - N) / N * pi
        y[i][j] = 2.0 * j / N * pi - pi / 2.0

Fa = np.zeros((N / 2 + 1, N))
Fb = np.zeros((N / 2 + 1, N))
T = np.zeros(N)

a_coef = np.random.normal(size=(L_max_polynom + 1, L_max_polynom + 1))
b_coef = np.random.normal(size=(L_max_polynom + 1, L_max_polynom + 1))

# a_coef = np.zeros((L_max_field + 1, L_max_field + 1))
# b_coef = np.zeros((L_max_field + 1, L_max_field + 1))

for m in xrange(0, L_max_field + 1):
    for l in xrange(0, m):
        a_coef[m][l] = 0.0
        b_coef[m][l] = 0.0
for l in xrange(0, L_max_field + 1):
    b_coef[0][l] = 0.0

a_coef[0][0] = 0.0
b_coef[0][0] = 0.0
a_coef[0][1] = 0.0
a_coef[1][1] = 0.0
b_coef[0][1] = 0.0
b_coef[1][1] = 0.0

C = np.zeros((L_max_field + 1))

for l in xrange(0, L_max_field + 1):

    C_sum = 0.0

    for m in xrange(0, l + 1):
        C_sum = C_sum + a_coef[m][l] * a_coef[m][l] + b_coef[m][l] * b_coef[m][l]

    C[l] = C_sum / (2.0 * l + 1.0)

sigma_0 = 0.0
for l in xrange(0, L_max_field + 1):
    sigma_0 += (2.0 * l + 1.0) * C[l]
sigma_0 = sqrt(sigma_0 / 4.0 / pi)

sigma_1 = 0.0
for l in xrange(0, L_max_field + 1):
    sigma_1 += l * (l + 1.0) * (2.0 * l + 1.0) * C[l]
sigma_1 = sqrt(sigma_1 / 4.0 * pi)

sigma_2 = 0.0
for l in xrange(0, L_max_field + 1):
    sigma_2 += (l + 2.0) * (l - 1.0) * l * (l + 1.0) * (2.0 * l + 1.0) * C[l]
sigma_2 = sqrt(sigma_2 / 4.0 * pi)

func1 = 0.0
func2 = 0.0

field = direct_f(N, P_, Fa, Fb, a_coef, b_coef, L_max_field)

field_x = direct_f(N, F_x, Fa, Fb, a_coef, b_coef, L_max_field, True)

field_y = direct_f(N, F_y, Fa, Fb, a_coef, b_coef, L_max_field)

field_xx = direct_f(N, F_xx_1 + F_xx_2, Fa, Fb, a_coef, b_coef, L_max_field)

field_yy = direct_f(N, F_yy, Fa, Fb, a_coef, b_coef, L_max_field)

field_xy = direct_f(N, F_xy, Fa, Fb, a_coef, b_coef, L_max_field, True)

a = 0.0
na = 0.0

for i in xrange(0, N):
    for j in xrange(1, N / 2):
        a += cos(y[i][j]) * field[i][j] * field[i][j]
        na += cos(y[i][j])

sigma_0_map = sqrt(a / na)

if sigma_0_map == 0:
    print 'There is no map!'

# field /= sigma_0_map

U = 2 * field_xy
Q = field_xx - field_yy
P = U*U + Q*Q

back_F = np.zeros((N / 2 + 1, N), dtype=complex)
back_F_a = np.zeros((N / 2 + 1, N))
back_F_b = np.zeros((N / 2 + 1, N))

back_a_coef, back_b_coef = inverse_f(N, field, P_, back_F, back_F_a, back_F_b, L_max_back)

my_map = cmbplot.flat(x, y, P)
cmbplot.polarization(N, my_map, U, Q, x, y)
cmbplot.show()
