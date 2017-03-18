from numpy import genfromtxt, zeros
from math import pi
from lib import cmbplot

N = 1024

file_map = genfromtxt('planck_2_dir/file_map_U.dat')

file_normal = zeros((N + 1, N / 2 + 1))

x = zeros((N + 1, N / 2 + 1))
y = zeros((N + 1, N / 2 + 1))

for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        x[i][j] = (2.0 * i - N) / N * pi
        y[i][j] = 2.0 * j / N * pi - pi / 2.0

for i in xrange(0, (N + 1) * (N/2 + 1)):
    file_normal[int(file_map[i][1]), int(file_map[i][2])] = file_map[i][0]

cmb_map = cmbplot.moll(x, y, file_normal)

cmbplot.level_plot(cmb_map, file_normal, x, y, 0.0)

cmbplot.show()