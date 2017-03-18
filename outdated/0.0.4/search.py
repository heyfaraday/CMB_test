from pylab import *
import numpy as np
import copy

x,y,z = genfromtxt('f1.dat').T
N = int(sqrt(x.shape[0]))

x = x.reshape(N, N)
y = y.reshape(N, N)
z = z.reshape(N, N)
mark1 = copy.deepcopy(z)
mark2 = copy.deepcopy(z)

for i in range (0, N-1):
    z_point = z[i] * z[i+1]
    for j in range (0, N):
        if (z_point[j] < 0) :
            mark1[i][j] = 1.0
        else:
            mark1[i][j] = 0.0
mark1[253] = mark1[253] * 0.0

z = z.T
x = x.T
y = y.T

for i in range (0, N-1):
    z_point = z[i] * z[i+1]
    for j in range (0, N):
        if (z_point[j] < 0) :
            mark2[i][j] = 1.0
        else:
            mark2[i][j] = 0.0

mark2[253] = mark2[253] * 0.0
mark2 = mark2.T

#print mark1[0][226] + mark1[0][227] + mark2[0][226] + mark2[1][226]

cell1 = copy.deepcopy(mark1)

for i in range (0, N-1):
    for j in range (0, N-1):
        cell1[i][j] = mark1[i][j] + mark1[i][j+1] + mark2[i][j] + mark2[i+1][j]

x,y,z = genfromtxt('f2.dat').T
N = int(sqrt(x.shape[0]))

x = x.reshape(N, N)
y = y.reshape(N, N)
z = z.reshape(N, N)
mark1 = copy.deepcopy(z)
mark2 = copy.deepcopy(z)
cell2 = copy.deepcopy(z)

for i in range (0, N-1):
    z_point = z[i] * z[i+1]
    for j in range (0, N):
        if (z_point[j] < 0) :
            mark1[i][j] = 1.0
        else:
            mark1[i][j] = 0.0
mark1[253] = mark1[253] * 0.0

z = z.T
x = x.T
y = y.T

for i in range (0, N-1):
    z_point = z[i] * z[i+1]
    for j in range (0, N):
        if (z_point[j] < 0) :
            mark2[i][j] = 1.0
        else:
            mark2[i][j] = 0.0

mark2[253] = mark2[253] * 0.0

mark2 = mark2.T

cell2 = copy.deepcopy(mark1)

for i in range (0, N-1):
    for j in range (0, N-1):
        cell2[i][j] = mark1[i][j] + mark1[i][j+1] + mark2[i][j] + mark2[i+1][j]

cell = cell1*cell2

answer_x = []
answer_y = []

for i in range (0, N-1):
    for j in range (0, N-1):
        if (cell[i][j] != 0): #& (cell[i][j] != 4):
            answer_x.append(float(i)/254)
            answer_y.append(float(j)/254)


x1,y1,z1 = genfromtxt('f1.dat').T
x1 = x.reshape(N,N)
y1 = y.reshape(N,N)
z1 = z.reshape(N,N)

x2,y2,z2 = genfromtxt('f2.dat').T
x2 = x.reshape(N,N)
y2 = y.reshape(N,N)
z2 = z.reshape(N,N)

pcolor(x,y,sqrt(z1*z1+z2*z2))

plot(answer_x, answer_y, 'kx', ms=20)

a1, b1 = genfromtxt('lev1.dat').T
plot(a1, b1, 'ko', ms=1)

a2, b2 = genfromtxt('lev2.dat').T
plot(a2, b2, 'ko', ms=1)

show()

file = open("points.dat", 'w')
for k in range(0, size(answer_x)-1):
    file.write('{}    {}\n'.format(int(254*answer_x[k]), int(254*answer_y[k])))

#-----------------------------------------

from numpy import roots
from math import atan, fabs

def cubic (Qx, Qy, Ux, Uy):
    
    a = Uy
    b = (Ux + 2*Qy)
    c = (2*Qx - Uy)
    d = -Ux
    
    det = -4*b*b*b*d + b*b*c*c -4*a*c*c*c + 18*a*b*c*d - 27*a*a*d*d
    
    print det
    
    if (det < 0):
        answer = 'c'
        #print answer
        
    if (det > 0):
        a = roots([a, b, c, d])
        a = a.real
        #print a
        a = atan(a[0])/(fabs(atan(a[0]))) + atan(a[1])/(fabs(atan(a[1]))) + atan(a[2])/(fabs(atan(a[2])))
        
        if (a == 3.0):
            answer = 'b'
        if (a == 1.0):
            answer = 'b'
        if (a == -1.0):
            answer = 'a'
        #print answer
        
    return answer


#-----------------------------------------

type = []

Ux, Uy, Uxx, Uyy, Uxy = genfromtxt('derivatives1.dat').T
Qx, Qy, Qxx, Qyy, Qxy = genfromtxt('derivatives2.dat').T

M = int(sqrt(Ux.shape[0]))

Ux = Ux.reshape(M, M)
Uy = Uy.reshape(M, M)

Qx = Qx.reshape(M, M)
Qy = Qy.reshape(M, M)

for k in range(0, size(answer_x)-1):
    type.append(cubic(Qx[answer_x[k]][answer_y[k]], Qy[answer_x[k]][answer_y[k]], Ux[answer_x[k]][answer_y[k]], Uy[answer_x[k]][answer_y[k]]))

print type 

file = open("points_types.dat", 'w')
for k in range(0, size(answer_x)-1):
    file.write('{}    {}    {}\n'.format(int(254*answer_x[k]), int(254*answer_y[k]), type[k]))
