import numpy as np
from math import sqrt, pi, sin, cos
from lib import cmbplot

L_max_field = 1
L_max_polynom = 1
N = 1024


def coef_1(in_l, in_m):
    if in_l != 0:
        return sqrt((in_l - in_m) * (2.0 * in_l + 1.0)
                    / ((in_l + in_m) * (2.0 * in_l - 1.0)))
    if in_l == 0:
        return 0.0


a_coef = np.random.normal(size=(L_max_field + 1, L_max_field + 1))
b_coef = np.random.normal(size=(L_max_field + 1, L_max_field + 1))

x = np.zeros((N + 1, N / 2 + 1))
y = np.zeros((N + 1, N / 2 + 1))

for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        x[i][j] = (2.0 * i - N) / N * pi
        y[i][j] = 2.0 * j / N * pi - pi / 2.0

Fa = np.zeros((N / 2 + 1, N))
Fb = np.zeros((N / 2 + 1, N))

func1 = 0.0
func2 = 0.0

P_ = np.zeros((L_max_polynom + 1, L_max_polynom + 1))

field = np.zeros((N + 1, N / 2 + 1))
for j in xrange(1, N / 2):

    theta = 2.0 * pi * j / N

    P_[0][0] = 1.0 / sqrt(4.0 * pi)

    for m in xrange(0, L_max_polynom):
        P_[m + 1][m + 1] = - P_[m][m] * sin(theta) * sqrt(2.0 * m + 3.0) / sqrt(2.0 * m + 2.0)

    for m in xrange(0, L_max_polynom):
        P_[m][m + 1] = P_[m][m] * cos(theta) * sqrt(2.0 * m + 3.0)

    for m in xrange(0, L_max_polynom - 1):
        for l in xrange(m + 2, L_max_polynom + 1):
            P_[m][l] = ((2.0 * l - 1.0) * sqrt((l - m) * (2.0 * l + 1.0)) / sqrt((l + m) * (2.0 * l - 1.0)) *
                            P_[m][l - 1] * cos(theta) - (l + m - 1.0) * sqrt((l - m) * (l - 1.0 - m) *
                            (2.0 * l + 1.0)) / sqrt((l + m) * (l - 1.0 + m) * (2.0 * l - 3.0)) * P_[m][l - 2]) / \
                            (l - m)

    for m in xrange(1, L_max_polynom + 1):
        for l in xrange(m, L_max_polynom + 1):
            P_[m][l] *= sqrt(2.0)

    for m in xrange(0, L_max_field + 1):
        for l in xrange(m, L_max_field + 1):
            func1 += a_coef[m][l] * P_[m][l]
            func2 += b_coef[m][l] * P_[m][l]

        Fa[j][m] = func1
        Fb[j][m] = func2

        func1 = 0.0
        func2 = 0.0

    T = np.real(np.fft.fft(Fa[j])) + np.imag(np.fft.fft(Fb[j]))

    field[0:N, j] = T[:]
    field[N][j] = field[0][j]

my_map = cmbplot.moll(x, y, field)
cmbplot.contour(my_map, field, x, y, 2)
cmbplot.show()
