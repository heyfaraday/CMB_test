from pylab import *
import numpy as np
import copy
from numpy.fft import ifft2
from copy import deepcopy

x,y,z = genfromtxt('f2.dat').T
N = int(sqrt(x.shape[0]))
x = x.reshape(N, N)
y = y.reshape(N, N)
z = z.reshape(N, N)

a1 = deepcopy(z)
b1 = deepcopy(z)

a1_x = deepcopy(z)
b1_x = deepcopy(z)
a1_y = deepcopy(z)
b1_y = deepcopy(z)

a1_xx = deepcopy(z)
b1_xx = deepcopy(z)
a1_yy = deepcopy(z)
b1_yy = deepcopy(z)
a1_xy = deepcopy(z)
b1_xy = deepcopy(z)

a1_xxx = deepcopy(z)
b1_xxx = deepcopy(z)
a1_yyy = deepcopy(z)
b1_yyy = deepcopy(z)
a1_xyy = deepcopy(z)
b1_xyy = deepcopy(z)
a1_xxy = deepcopy(z)
b1_xxy = deepcopy(z)


a2 = deepcopy(z)
b2 = deepcopy(z)

a2_x = deepcopy(z)
b2_x = deepcopy(z)
a2_y = deepcopy(z)
b2_y = deepcopy(z)

a2_xx = deepcopy(z)
b2_xx = deepcopy(z)
a2_yy = deepcopy(z)
b2_yy = deepcopy(z)
a2_xy = deepcopy(z)
b2_xy = deepcopy(z)

a2_xxx = deepcopy(z)
b2_xxx = deepcopy(z)
a2_yyy = deepcopy(z)
b2_yyy = deepcopy(z)
a2_xyy = deepcopy(z)
b2_xyy = deepcopy(z)
a2_xxy = deepcopy(z)
b2_xxy = deepcopy(z)


modes = 7;

mu, sigma = 0.0, 1.0 

def P(i_f, j_f):
    if (i_f > modes) | (i_f < -modes): #
        return 0.0 #
    if (j_f < 0) | (j_f > modes): #
        return 0.0 #

    if (j_f > N/2):
        return 0.0
    if (i_f**2 + j_f**2 > modes**2):
        return 0
    if (j_f == 0) & (i_f == 0):
        return 0
    return 1
    #return exp(-(i_f**2 + j_f**2)/(modes))

for i in xrange(0, N):
    for j in xrange(0, N):
        i_f = float(i-N/2)
        j_f = float(j)
        a = np.random.normal(mu, sigma)
        b = np.random.normal(mu, sigma)

        a1[int(i)][int(j)] = float(P(i_f, j_f)*(-2*i_f*j_f)*a)
        b1[int(i)][int(j)] = float(P(i_f, j_f)*(-2*i_f*j_f)*b)

        a1_x[int(i)][int(j)] = i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*a)
        b1_x[int(i)][int(j)] = i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*b)
        a1_y[int(i)][int(j)] = j_f*float(P(i_f, j_f)*(-2*i_f*j_f)*a)
        b1_y[int(i)][int(j)] = j_f*float(P(i_f, j_f)*(-2*i_f*j_f)*b)

        a1_xx[int(i)][int(j)] = i_f*i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*a)
        b1_xx[int(i)][int(j)] = i_f*i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*b)
        a1_yy[int(i)][int(j)] = j_f*j_f*float(P(i_f, j_f)*(-2*i_f*j_f)*a)
        b1_yy[int(i)][int(j)] = j_f*j_f*float(P(i_f, j_f)*(-2*i_f*j_f)*b)
        a1_xy[int(i)][int(j)] = j_f*i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*a)
        b1_xy[int(i)][int(j)] = j_f*i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*b)

        a1_xxx[int(i)][int(j)] = i_f*i_f*i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*a)
        b1_xxx[int(i)][int(j)] = i_f*i_f*i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*b)
        a1_yyy[int(i)][int(j)] = j_f*j_f*j_f*float(P(i_f, j_f)*(-2*i_f*j_f)*a)
        b1_yyy[int(i)][int(j)] = j_f*j_f*j_f*float(P(i_f, j_f)*(-2*i_f*j_f)*b)
        a1_xyy[int(i)][int(j)] = j_f*j_f*i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*a)
        b1_xyy[int(i)][int(j)] = j_f*j_f*i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*b)
        a1_xxy[int(i)][int(j)] = j_f*i_f*i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*a)
        b1_xxy[int(i)][int(j)] = j_f*i_f*i_f*float(P(i_f, j_f)*(-2*i_f*j_f)*b)

for i in xrange(0, N):
    for j in xrange(0, N):
        i_f = float(i-N/2)
        j_f = float(j)
        a = np.random.normal(mu, sigma)
        b = np.random.normal(mu, sigma)

        a2[int(i)][int(j)] = float(P(i_f, j_f)*(-i_f**2+j_f**2)*a)
        b2[int(i)][int(j)] = float(P(i_f, j_f)*(-i_f**2+j_f**2)*b)

        a2_x[int(i)][int(j)] = i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*a)
        b2_x[int(i)][int(j)] = i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*b)
        a2_y[int(i)][int(j)] = j_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*a)
        b2_y[int(i)][int(j)] = j_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*b)

        a2_xx[int(i)][int(j)] = i_f*i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*a)
        b2_xx[int(i)][int(j)] = i_f*i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*b)
        a2_yy[int(i)][int(j)] = j_f*j_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*a)
        b2_yy[int(i)][int(j)] = j_f*j_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*b)
        a2_xy[int(i)][int(j)] = j_f*i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*a)
        b2_xy[int(i)][int(j)] = j_f*i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*b)

        a2_xxx[int(i)][int(j)] = i_f*i_f*i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*a)
        b2_xxx[int(i)][int(j)] = i_f*i_f*i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*b)
        a2_yyy[int(i)][int(j)] = j_f*j_f*j_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*a)
        b2_yyy[int(i)][int(j)] = j_f*j_f*j_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*b)
        a2_xyy[int(i)][int(j)] = j_f*j_f*i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*a)
        b2_xyy[int(i)][int(j)] = j_f*j_f*i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*b)
        a2_xxy[int(i)][int(j)] = j_f*i_f*i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*a)
        b2_xxy[int(i)][int(j)] = j_f*i_f*i_f*float(P(i_f, j_f)*(-i_f**2+j_f**2)*b)

U = np.real(fft2(a1)) + np.imag(fft2(b1))

U_x = -np.imag(fft2(a1_x)) + np.real(fft2(b1_x))
U_y = -np.imag(fft2(a1_y)) + np.real(fft2(b1_y))

U_xx = -np.real(fft2(a1_xx)) - np.imag(fft2(b1_xx))
U_yy = -np.real(fft2(a1_yy)) - np.imag(fft2(b1_yy))
U_xy = -np.real(fft2(a1_xy)) - np.imag(fft2(b1_xy))

U_xxx = np.imag(fft2(a1_xxx)) - np.real(fft2(b1_xxx))
U_yyy = np.imag(fft2(a1_yyy)) - np.real(fft2(b1_yyy))
U_xyy = np.imag(fft2(a1_xyy)) - np.real(fft2(b1_xyy))
U_xxy = np.imag(fft2(a1_xxy)) - np.real(fft2(b1_xxy))


Q = np.real(fft2(a2)) + np.imag(fft2(b2))

Q_x = -np.imag(fft2(a2_x)) + np.real(fft2(b2_x))
Q_y = -np.imag(fft2(a2_y)) + np.real(fft2(b2_y))

Q_xx = -np.real(fft2(a2_xx)) - np.imag(fft2(b2_xx))
Q_yy = -np.real(fft2(a2_yy)) - np.imag(fft2(b2_yy))
Q_xy = -np.real(fft2(a2_xy)) - np.imag(fft2(b2_xy))

Q_xxx = np.imag(fft2(a2_xxx)) - np.real(fft2(b2_xxx))
Q_yyy = np.imag(fft2(a2_yyy)) - np.real(fft2(b2_yyy))
Q_xyy = np.imag(fft2(a2_xyy)) - np.real(fft2(b2_xyy))
Q_xxy = np.imag(fft2(a2_xxy)) - np.real(fft2(b2_xxy))

for i in xrange(0, N):
    for j in xrange(0, N):
        if (i % 2 == 1):
            U[i][j] = -U[i][j]

            U_x[i][j] = -U_x[i][j]
            U_y[i][j] = -U_y[i][j]

            U_xx[i][j] = -U_xx[i][j]
            U_yy[i][j] = -U_yy[i][j]
            U_xy[i][j] = -U_xy[i][j]

            U_xxx[i][j] = -U_xxx[i][j]
            U_yyy[i][j] = -U_yyy[i][j]
            U_xyy[i][j] = -U_xyy[i][j]
            U_xxy[i][j] = -U_xxy[i][j]

            Q[i][j] = -Q[i][j]

            Q_x[i][j] = -Q_x[i][j]
            Q_y[i][j] = -Q_y[i][j]

            Q_xx[i][j] = -Q_xx[i][j]
            Q_yy[i][j] = -Q_yy[i][j]
            Q_xy[i][j] = -Q_xy[i][j]

            Q_xxx[i][j] = -Q_xxx[i][j]
            Q_yyy[i][j] = -Q_yyy[i][j]
            Q_xyy[i][j] = -Q_xyy[i][j]
            Q_xxy[i][j] = -Q_xxy[i][j]

pcolormesh(x,y,U_xx)
show()

pcolormesh(x,y,U_yy)
show()

pcolormesh(x,y,np.sqrt(np.sqrt(Q**2 + U**2)))
show()

file = open('U.dat', 'w')
for j1 in range(0, N):
    for j2 in range (0, N):
        x = float(j1)/N
        y = float(j2)/N
        file.write(repr(U[j1][j2]) + '  ' + repr(x) + '   ' + repr(y) + '\n')

file = open('derivativesU.dat', 'w')
for j1 in range(0, N):
    for j2 in range (0, N):
        file.write(repr(U_x[j1][j2]) + '  ' + repr(U_y[j1][j2]) + '   ' + repr(U_xx[j1][j2]) + '   ' + repr(U_yy[j1][j2]) + '   ' + repr(U_xy[j1][j2]) + '\n')

file = open('Q.dat', 'w')
for j1 in range(0, N):
    for j2 in range (0, N):
        x = float(j1)/N
        y = float(j2)/N
        file.write(repr(Q[j1][j2]) + '  ' + repr(x) + '   ' + repr(y) + '\n')

file = open('derivativesQ.dat', 'w')
for j1 in range(0, N):
    for j2 in range (0, N):
        file.write(repr(Q_x[j1][j2]) + '  ' + repr(Q_y[j1][j2]) + '   ' + repr(Q_xx[j1][j2]) + '   ' + repr(Q_yy[j1][j2]) + '   ' + repr(Q_xy[j1][j2]) + '\n')

file = open('skfU.dat', 'w')
for j1 in range(0, N):
    for j2 in range (0, N):
        dens = U_xx[j1][j2] + U_yy[j1][j2]
        dfs = U_x[j1][j2]*U_y[j1][j2]*(U_yy[j1][j2] - U_xx[j1][j2]) + (U_x[j1][j2]**2-U_y[j1][j2]**2)*U_xy[j1][j2]
        dfss = U_y[j1][j2]*U_y[j1][j2]*(U_xx[j1][j2]**2+U_x[j1][j2]*U_xxx[j1][j2]+U_xy[j1][j2]**2+U_y[j1][j2]*U_xxy[j1][j2]) + U_x[j1][j2]*U_x[j1][j2]*(U_yy[j1][j2]**2+U_y[j1][j2]*U_yyy[j1][j2]+U_xy[j1][j2]**2+U_x[j1][j2]*U_xyy[j1][j2]) - 2*U_x[j1][j2]*U_y[j1][j2]*(U_xy[j1][j2]*(U_xx[j1][j2]+U_yy[j1][j2]) + U_x[j1][j2]*U_xxy[j1][j2]+U_y[j1][j2]*U_xyy[j1][j2])
        file.write(repr(dens) + '  ' + repr(dfs) + '   ' + repr(dfss) + '\n')


#file = open('P.dat', 'w')
#for j1 in range(0, N):
#    for j2 in range (0, N):
#        x = float(j1)/N
#        y = float(j2)/N
#        file.write(repr(sqrt(Q[j1][j2]**2 + U[j1][j2]**2)) + '  ' + repr(x) + '   ' + repr(y) + '\n')

#file = open('derivativesP.dat', 'w')
#for j1 in range(0, N):
#    for j2 in range (0, N):
#        x = float(j1)/N
#        y = float(j2)/N
#        P_x = 
#        file.write(repr(Q_x[j1][j2]) + '  ' + repr(Q_y[j1][j2]) + '   ' + repr(Q_xx[j1][j2]) + repr(Q_yy[j1][j2]) + repr(Q_xy[j1][j2]) + '\n')
