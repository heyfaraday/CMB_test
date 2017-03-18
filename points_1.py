from numpy import genfromtxt, zeros
from math import pi
from lib import minkowski

N = 1024

file_map_Q = genfromtxt('planck_2_dir/file_map_Q.dat')
file_map_Qx = genfromtxt('planck_2_dir/file_map_Qx.dat')
file_map_Qy = genfromtxt('planck_2_dir/file_map_Qy.dat')

file_map_U = genfromtxt('planck_2_dir/file_map_U.dat')
file_map_Ux = genfromtxt('planck_2_dir/file_map_Ux.dat')
file_map_Uy = genfromtxt('planck_2_dir/file_map_Uy.dat')

file_map_points = open('planck_2_dir/file_map_points.dat', 'w')

file_normal_Q = zeros((N + 1, N / 2 + 1))
file_normal_Qx = zeros((N + 1, N / 2 + 1))
file_normal_Qy = zeros((N + 1, N / 2 + 1))
file_normal_U = zeros((N + 1, N / 2 + 1))
file_normal_Ux = zeros((N + 1, N / 2 + 1))
file_normal_Uy = zeros((N + 1, N / 2 + 1))


x = zeros((N + 1, N / 2 + 1))
y = zeros((N + 1, N / 2 + 1))

for i in xrange(0, N + 1):
    for j in xrange(0, N / 2 + 1):
        x[i][j] = (2.0 * i - N) / N * pi
        y[i][j] = 2.0 * j / N * pi - pi / 2.0

for i in xrange(0, (N + 1) * (N/2 + 1)):
    file_normal_Q[int(file_map_Q[i][1]), int(file_map_Q[i][2])] = file_map_Q[i][0]
    file_normal_Qx[int(file_map_Qx[i][1]), int(file_map_Qx[i][2])] = file_map_Qx[i][0]
    file_normal_Qy[int(file_map_Qy[i][1]), int(file_map_Qy[i][2])] = file_map_Qy[i][0]
    file_normal_U[int(file_map_U[i][1]), int(file_map_U[i][2])] = file_map_U[i][0]
    file_normal_Ux[int(file_map_Ux[i][1]), int(file_map_Ux[i][2])] = file_map_Ux[i][0]
    file_normal_Uy[int(file_map_Uy[i][1]), int(file_map_Uy[i][2])] = file_map_Uy[i][0]


# cmb_map = cmbplot.moll(x, y, file_normal_Q*file_normal_Q + file_normal_U*file_normal_U)

minkowski.singular_points(x, y, file_normal_Q, file_normal_U, file_normal_Qx, file_normal_Qy, file_normal_Ux,
                          file_normal_Uy, my_file=file_map_points)
