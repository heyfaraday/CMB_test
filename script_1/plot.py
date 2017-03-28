from numpy import genfromtxt, zeros, size, sqrt
from math import pi
from lib import cmbplot, minkowski

file_map_Q_512 = genfromtxt('data/map_Q_64_128.dat')
file_map_U_512 = genfromtxt('data/map_U_64_128.dat')
file_map_points_512 = genfromtxt('data/file_map_points_64_128.dat')

N = 128
n_points_512 = int(size(file_map_points_512) / 9.0)

file_normal_Q_512 = zeros((N + 1, N / 2 + 1))
file_normal_U_512 = zeros((N + 1, N / 2 + 1))

x = zeros((N + 1, N / 2 + 1))
y = zeros((N + 1, N / 2 + 1))

for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        x[i][j] = (2.0 * i - N) / N * pi
        y[i][j] = 2.0 * j / N * pi - pi / 2.0

for i in xrange(0, (N + 1) * (N / 2 + 1)):
    file_normal_Q_512[int(file_map_Q_512[i][1]), int(file_map_Q_512[i][2])] = file_map_Q_512[i][0]
    file_normal_U_512[int(file_map_U_512[i][1]), int(file_map_U_512[i][2])] = file_map_U_512[i][0]

P = sqrt(file_normal_Q_512 * file_normal_Q_512 + file_normal_U_512 * file_normal_U_512)

cmb_map = cmbplot.flat(x, y, P)
cmbplot.contour(cmb_map, file_normal_Q_512, x, y, 1)
cmbplot.contour(cmb_map, file_normal_U_512, x, y, 1)

minkowski.points_comparison_single(file_map_points_512, N, my_cmbmap=cmb_map, number_plot=3000, pix=True)

cmbplot.show()
