from numpy import genfromtxt, zeros, size, sqrt
from math import pi
from lib import cmbplot, minkowski

file_map_Q_512 = genfromtxt('planck_2_dir/file_map_Q_64_4096.dat')
file_map_U_512 = genfromtxt('planck_2_dir/file_map_U_64_4096.dat')

file_map_points_512 = genfromtxt('planck_2_dir/file_map_points_64_4096.dat')

file_map_Q_1024 = genfromtxt('planck_2_dir/file_map_Q_64_4096.dat')
file_map_U_1024 = genfromtxt('planck_2_dir/file_map_U_64_4096.dat')

file_map_points_1024 = genfromtxt('planck_2_dir/file_map_points_64_4096.dat')

N = 2048
n_points_512 = int(size(file_map_points_512) / 9.0)
n_points_1024 = int(size(file_map_points_1024) / 9.0)

file_normal_Q_512 = zeros((N + 1, N / 2 + 1))
file_normal_U_512 = zeros((N + 1, N / 2 + 1))

file_normal_Q_1024 = zeros((N + 1, N / 2 + 1))
file_normal_U_1024 = zeros((N + 1, N / 2 + 1))

x = zeros((N + 1, N / 2 + 1))
y = zeros((N + 1, N / 2 + 1))

for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        x[i][j] = (2.0 * i - N) / N * pi
        y[i][j] = 2.0 * j / N * pi - pi / 2.0

for i in xrange(0, (N + 1) * (N/2 + 1)):
    file_normal_Q_512[int(file_map_Q_512[i][1]), int(file_map_Q_512[i][2])] = file_map_Q_512[i][0]
    file_normal_U_512[int(file_map_U_512[i][1]), int(file_map_U_512[i][2])] = file_map_U_512[i][0]
    file_normal_Q_1024[int(file_map_Q_1024[i][1]), int(file_map_Q_1024[i][2])] = file_map_Q_1024[i][0]
    file_normal_U_1024[int(file_map_U_1024[i][1]), int(file_map_U_1024[i][2])] = file_map_U_1024[i][0]

cmb_map = cmbplot.flat(x, y, sqrt(file_normal_Q_1024 * file_normal_Q_1024 + file_normal_U_1024 * file_normal_U_1024))

# cmbplot.level_plot(cmb_map, file_normal_Q_1024, x, y, 0.0)
# cmbplot.level_plot(cmb_map, file_normal_U_1024, x, y, 0.0)

# minkowski.null_points(x, y, sqrt(file_normal_Q_1024 * file_normal_Q_1024 + file_normal_U_1024 * file_normal_U_1024),
#                       file_normal_U_1024, file_normal_Q_1024, )

minkowski.points_comparison(file_map_points_512, file_map_points_512, N, my_cmbmap=cmb_map, number_plot=3000)

cmbplot.show()
