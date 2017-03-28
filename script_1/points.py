from numpy import genfromtxt, zeros, real, imag
from math import pi
from lib import minkowski, healpy_format

N = 128
L_max_field = 64

file_map_Q = genfromtxt('data/map_Q_64_128.dat')
file_map_Qx = genfromtxt('data/map_Qx_64_128.dat')
file_map_Qy = genfromtxt('data/map_Qy_64_128.dat')

file_map_U = genfromtxt('data/map_U_64_128.dat')
file_map_Ux = genfromtxt('data/map_Ux_64_128.dat')
file_map_Uy = genfromtxt('data/map_Uy_64_128.dat')

file_map_points = open('data/file_map_points_64_128.dat', 'w')

file_normal_Q = zeros((N + 1, N / 2 + 1))
file_normal_Qx = zeros((N + 1, N / 2 + 1))
file_normal_Qy = zeros((N + 1, N / 2 + 1))
file_normal_U = zeros((N + 1, N / 2 + 1))
file_normal_Ux = zeros((N + 1, N / 2 + 1))
file_normal_Uy = zeros((N + 1, N / 2 + 1))

alm_file_Q = healpy_format.healpy_file('data/alm_Q_norm_64.dat', L_max_field)
a_coef_Q = real(alm_file_Q)
b_coef_Q = imag(alm_file_Q)

alm_file_U = healpy_format.healpy_file('data/alm_U_norm_64.dat', L_max_field)
a_coef_U = real(alm_file_U)
b_coef_U = imag(alm_file_U)

a_coef_Q[0][0] = 0.0
b_coef_Q[0][0] = 0.0
a_coef_Q[0][1] = 0.0
a_coef_Q[1][1] = 0.0
b_coef_Q[0][1] = 0.0
b_coef_Q[1][1] = 0.0

x = zeros((N + 1, N / 2 + 1))
y = zeros((N + 1, N / 2 + 1))

for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        x[i][j] = (2.0 * i - N) / N * pi
        y[i][j] = 2.0 * j / N * pi - pi / 2.0

for i in xrange(0, (N + 1) * (N / 2 + 1)):
    file_normal_Q[int(file_map_Q[i][1]), int(file_map_Q[i][2])] = file_map_Q[i][0]
    file_normal_Qx[int(file_map_Qx[i][1]), int(file_map_Qx[i][2])] = file_map_Qx[i][0]
    file_normal_Qy[int(file_map_Qy[i][1]), int(file_map_Qy[i][2])] = file_map_Qy[i][0]
    file_normal_U[int(file_map_U[i][1]), int(file_map_U[i][2])] = file_map_U[i][0]
    file_normal_Ux[int(file_map_Ux[i][1]), int(file_map_Ux[i][2])] = file_map_Ux[i][0]
    file_normal_Uy[int(file_map_Uy[i][1]), int(file_map_Uy[i][2])] = file_map_Uy[i][0]

n_saddle, n_beak, n_comet = minkowski.singular_points(x, y, file_normal_Q, file_normal_U, file_normal_Qx,
                            file_normal_Qy, file_normal_Ux, file_normal_Uy, a_coef_Q, b_coef_Q,  a_coef_U,
                            b_coef_U, L_max_field, my_file=file_map_points, print_num=True)

print 'saddles:', n_saddle
print 'beaks:', n_beak
print 'comets', n_comet
