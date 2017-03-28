from math import sqrt, pi, sin, cos
import numpy as np
from lib.fourier import direct_f, direct_f_int
from lib.minkowski import length, area, type_points, null_points
from lib import cmbplot
from lib import healpy_format

from lib.fourier import direct_point_int, inverse_f_int, inverse_f_spin_int

from numpy.fft import fft
from numpy import zeros

L_max_field = 2047
L_max_polynom = 2047
L_max_back = 2047
N = 2048

# file_map_I = open('data/map_I_2048_2048.dat', 'w')
# file_map_Ix = open('data/map_Ix_2048_2048.dat', 'w')
# file_map_Iy = open('data/map_Iy_2048_2048.dat', 'w')

file_map_Q = open('data/map_Q_2048_2048.dat', 'w')
file_map_Qx = open('data/map_Qx_2048_2048.dat', 'w')
file_map_Qy = open('data/map_Qy_2048_2048.dat', 'w')

file_map_U = open('data/map_U_2048_2048.dat', 'w')
file_map_Ux = open('data/map_Ux_2048_2048.dat', 'w')
file_map_Uy = open('data/map_Uy_2048_2048.dat', 'w')

x = np.zeros((N + 1, N / 2 + 1))
y = np.zeros((N + 1, N / 2 + 1))

for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        x[i][j] = (2.0 * i - N) / N * pi
        y[i][j] = 2.0 * j / N * pi - pi / 2.0

# alm_file_I = healpy_format.healpy_file('data/alm_I_norm_2048.dat', L_max_field)
alm_file_Q = healpy_format.healpy_file('data/alm_Q_norm_2048.dat', L_max_field)
alm_file_U = healpy_format.healpy_file('data/alm_U_norm_2048.dat', L_max_field)

# a_coef_I = np.real(alm_file_I)
# b_coef_I = np.imag(alm_file_I)
a_coef_Q = np.real(alm_file_Q)
b_coef_Q = np.imag(alm_file_Q)
a_coef_U = np.real(alm_file_U)
b_coef_U = np.imag(alm_file_U)

# a_coef_I[0][0] = 0.0
# b_coef_I[0][0] = 0.0
# a_coef_I[0][1] = 0.0
# a_coef_I[1][1] = 0.0
# b_coef_I[0][1] = 0.0
# b_coef_I[1][1] = 0.0

a_coef_Q[0][0] = 0.0
b_coef_Q[0][0] = 0.0
a_coef_Q[0][1] = 0.0
a_coef_Q[1][1] = 0.0
b_coef_Q[0][1] = 0.0
b_coef_Q[1][1] = 0.0

a_coef_U[0][0] = 0.0
b_coef_U[0][0] = 0.0
a_coef_U[0][1] = 0.0
a_coef_U[1][1] = 0.0
b_coef_U[0][1] = 0.0
b_coef_U[1][1] = 0.0

# print 'I'
# I = direct_f_int(N, a_coef_I, b_coef_I, L_max_field, sign='0')
# print 'Ix'
# Ix = direct_f_int(N, a_coef_I, b_coef_I, L_max_field, sign='x', diff=True)
# print 'Iy'
# Iy = direct_f_int(N, a_coef_I, b_coefI, L_max_field, sign='y')

print 'Q'
Q = direct_f_int(N, a_coef_Q, b_coef_Q, L_max_field)
print 'Qx'
Qx = direct_f_int(N, a_coef_Q, b_coef_Q, L_max_field, sign=1, diff=True)
print 'Qy'
Qy = direct_f_int(N, a_coef_Q, b_coef_Q, L_max_field, sign=2)

print 'U'
U = direct_f_int(N, a_coef_U, b_coef_U, L_max_field)
print 'Ux'
Ux = direct_f_int(N, a_coef_U, b_coef_U, L_max_field, sign=1, diff=True)
print 'Uy'
Uy = direct_f_int(N, a_coef_U, b_coef_U, L_max_field, sign=2)


# print 'I_output'
# for i in xrange(0, N + 1):
#     for j in xrange(0, N / 2 + 1):
#         file_map_I.write(repr(I[i][j]) + '    ' + repr(i) + '    ' + repr(j) + '\n')
# for i in xrange(0, N + 1):
#     for j in xrange(0, N / 2 + 1):
#         file_map_Ix.write(repr(Ix[i][j]) + '    ' + repr(i) + '    ' + repr(j) + '\n')
# for i in xrange(0, N + 1):
#     for j in xrange(0, N / 2 + 1):
#         file_map_Iy.write(repr(Iy[i][j]) + '    ' + repr(i) + '    ' + repr(j) + '\n')

print 'Q_output'
for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        file_map_Q.write(repr(Q[i][j]) + '    ' + repr(i) + '    ' + repr(j) + '\n')
for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        file_map_Qx.write(repr(Qx[i][j]) + '    ' + repr(i) + '    ' + repr(j) + '\n')
for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        file_map_Qy.write(repr(Qy[i][j]) + '    ' + repr(i) + '    ' + repr(j) + '\n')

print 'U_output'
for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        file_map_U.write(repr(U[i][j]) + '    ' + repr(i) + '    ' + repr(j) + '\n')
for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        file_map_Ux.write(repr(Ux[i][j]) + '    ' + repr(i) + '    ' + repr(j) + '\n')
for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        file_map_Uy.write(repr(Uy[i][j]) + '    ' + repr(i) + '    ' + repr(j) + '\n')

# back Q
#
# Q_plus_U = Q.astype(complex)
# Q_plus_U.imag = U
#
# Q_minus_U = Q.astype(complex)
# Q_minus_U.imag = -U
#
# a_plus, b_plus = inverse_f_spin_int(N, Q_plus_U, 64)
# a_minus, b_minus = inverse_f_spin_int(N, Q_minus_U, 64, sign=-2)
#
# a_e = - (a_plus + a_minus) / 2.0
# a_b = - (a_plus + a_minus) / 2.0j
#
# b_e = - (b_plus + b_minus) / 2.0
# b_b = - (b_plus + b_minus) / 2.0j
#
# print a_plus[20][20]
# print a_minus[20][20]
