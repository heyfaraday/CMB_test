from math import sqrt, pi, sin, cos
import numpy as np
from lib.fourier import direct_f, direct_f_int
from lib.minkowski_1 import length, area, type_points, null_points
from lib import cmbplot
from lib import healpy_format

L_max_field = 512
L_max_polynom = 512
L_max_back = 512
N = 1024

file_map = open('planck_2_dir/file_map.dat', 'w')

genus = np.zeros(41)
Nmax = np.zeros(41)
Nmin = np.zeros(41)
Nsad = np.zeros(41)


def coef_1(in_l, in_m):
    if in_l != 0:
        return sqrt((in_l - in_m) * (2.0 * in_l + 1.0)
                    / ((in_l + in_m) * (2.0 * in_l - 1.0)))
    if in_l == 0:
        return 0.0

x = np.zeros((N + 1, N / 2 + 1))
y = np.zeros((N + 1, N / 2 + 1))

for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        x[i][j] = (2.0 * i - N) / N * pi
        y[i][j] = 2.0 * j / N * pi - pi / 2.0

Fa = np.zeros((N / 2 + 1, N))
Fb = np.zeros((N / 2 + 1, N))
T = np.zeros(N)

a_from_file = healpy_format.healpy_file('planck_2_dir/planck_2_I_norm.dat', L_max_field)

a_coef = np.real(a_from_file)
b_coef = np.imag(a_from_file)

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

field = direct_f_int(N, Fa, Fb, a_coef, b_coef, L_max_field, sign='0')
field_x = direct_f_int(N, Fa, Fb, a_coef, b_coef, L_max_field, sign='x', diff=True)
field_y = direct_f_int(N, Fa, Fb, a_coef, b_coef, L_max_field, sign='y')

a = 0.0
na = 0.0

for i in xrange(0, N):
    for j in xrange(1, N / 2):
        a += cos(y[i][j]) * field[i][j] * field[i][j]
        na += cos(y[i][j])

sigma_0_map = sqrt(a / na)

if sigma_0_map == 0:
    print 'There is no map!'

field /= sigma_0_map

for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        file_map.write(repr(field[i][j]) + '    ' + repr(i) + '    ' + repr(j) + '\n')
