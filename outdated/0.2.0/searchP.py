# Code works stable with environment.yml 

from pylab import *
import numpy as np
import copy

z,x,y = genfromtxt('U.dat').T
N = int(sqrt(x.shape[0]))

x = x.reshape(N, N)
y = y.reshape(N, N)
z = z.reshape(N, N)
mark_11 = copy.deepcopy(z)
mark_12 = copy.deepcopy(z)

for i in range (0, N-1):
    z_point = z[i] * z[i+1]
    for j in range (0, N):
        if (z_point[j] < 0) :
            mark_11[i][j] = 1.0
        else:
            mark_11[i][j] = 0.0
mark_11[253] = mark_11[253] * 0.0

z = z.T
x = x.T
y = y.T

for i in range (0, N-1):
    z_point = z[i] * z[i+1]
    for j in range (0, N):
        if (z_point[j] < 0) :
            mark_12[i][j] = 1.0
        else:
            mark_12[i][j] = 0.0

mark_12[253] = mark_12[253] * 0.0
mark_12 = mark_12.T

#print mark1[0][226] + mark1[0][227] + mark2[0][226] + mark2[1][226]

cell1 = copy.deepcopy(mark_11)

for i in range (0, N-1):
    for j in range (0, N-1):
        cell1[i][j] = mark_11[i][j] + mark_11[i][j+1] + mark_12[i][j] + mark_12[i+1][j]

z,x,y = genfromtxt('Q.dat').T
N = int(sqrt(x.shape[0]))

x = x.reshape(N, N)
y = y.reshape(N, N)
z = z.reshape(N, N)
mark_21 = copy.deepcopy(z)
mark_22 = copy.deepcopy(z)
cell2 = copy.deepcopy(z)

for i in range (0, N-1):
    z_point = z[i] * z[i+1]
    for j in range (0, N):
        if (z_point[j] < 0) :
            mark_21[i][j] = 1.0
        else:
            mark_21[i][j] = 0.0
mark_21[253] = mark_21[253] * 0.0

z = z.T
x = x.T
y = y.T

for i in range (0, N-1):
    z_point = z[i] * z[i+1]
    for j in range (0, N):
        if (z_point[j] < 0) :
            mark_22[i][j] = 1.0
        else:
            mark_22[i][j] = 0.0

mark_22[253] = mark_22[253] * 0.0

mark_22 = mark_22.T

cell2 = copy.deepcopy(mark_21)

for i in range (0, N-1):
    for j in range (0, N-1):
        cell2[i][j] = mark_21[i][j] + mark_21[i][j+1] + mark_22[i][j] + mark_22[i+1][j]

cell = cell1*cell2

answer_x = []
answer_y = []
answer_Ux = []
answer_Uy = []
answer_Qx = []
answer_Qy = []
 
z1,x1,y1 = genfromtxt('U.dat').T
x1 = x1.reshape(N,N)
y1 = y1.reshape(N,N)
z1 = z1.reshape(N,N)

z2,x2,y2 = genfromtxt('Q.dat').T
x2 = x2.reshape(N,N)
y2 = y2.reshape(N,N)
z2 = z2.reshape(N,N)

z1 = z1.T
z2 = z2.T


Ux, Uy, Uxx, Uyy, Uxy = genfromtxt('derivativesU.dat').T
Qx, Qy, Qxx, Qyy, Qxy = genfromtxt('derivativesQ.dat').T

M = int(sqrt(Ux.shape[0]))

Ux = Ux.reshape(M, M)
Uy = Uy.reshape(M, M)

Qx = Qx.reshape(M, M)
Qy = Qy.reshape(M, M)


for i in range (0, N-1):
    for j in range (0, N-1):
        if (cell[i][j] != 0): #& (cell[i][j] != 16):

            Ux1 = 0.0
            Uy1 = 0.0
            Ux2 = 0.0
            Uy2 = 0.0
            Qx1 = 0.0
            Qy1 = 0.0
            Qx2 = 0.0
            Qy2 = 0.0

            Uxx1 = 0.0
            Uyy1 = 0.0
            Uxx2 = 0.0
            Uyy2 = 0.0
            Qxx1 = 0.0
            Qyy1 = 0.0
            Qxx2 = 0.0
            Qyx2 = 0.0    

            if (mark_11[i][j] != 0.0):
                if (mark_11[i][j+1] != 0.0):
                    Ux1 = 0.0
                    Uy1 = fabs(float(z1[i][j]))/(fabs(float(z1[i][j]))+fabs(float(z1[i+1][j])))
                    Uxx1 = (Ux[i][j] - Ux[i+1][j])*Uy1 + Ux[i+1][j]
                    Uyy1 = (Uy[i][j] - Uy[i+1][j])*Uy1 + Uy[i+1][j]

                    Ux2 = 1.0
                    Uy2 = fabs(float(z1[i][j+1]))/(fabs(float(z1[i][j+1]))+fabs(float(z1[i+1][j+1])))
                    Uxx2 = (Ux[i][j+1] - Ux[i+1][j+1])*Uy2 + Ux[i+1][j+1]
                    Uyy2 = (Uy[i][j+1] - Uy[i+1][j+1])*Uy2 + Uy[i+1][j+1]

                if (mark_12[i][j] != 0.0):
                    Ux1 = 0.0
                    Uy1 = fabs(float(z1[i][j]))/(fabs(float(z1[i][j]))+fabs(float(z1[i+1][j])))
                    Uxx1 = (Ux[i][j] - Ux[i+1][j])*Uy1 + Ux[i+1][j]
                    Uyy1 = (Uy[i][j] - Uy[i+1][j])*Uy1 + Uy[i+1][j]

                    Uy2 = 0.0
                    Ux2 = fabs(float(z1[i][j]))/(fabs(float(z1[i][j+1]))+fabs(float(z1[i][j])))
                    Uxx2 = (Ux[i][j+1] - Ux[i][j])*Ux2 + Ux[i][j]
                    Uyy2 = (Uy[i][j+1] - Uy[i][j])*Ux2 + Uy[i][j]

                if (mark_12[i+1][j] != 0.0):
                    Ux1 = 0.0
                    Uy1 = fabs(float(z1[i][j]))/(fabs(float(z1[i][j]))+fabs(float(z1[i+1][j])))
                    Uxx1 = (Ux[i][j] - Ux[i+1][j])*Uy1 + Ux[i+1][j]
                    Uyy1 = (Uy[i][j] - Uy[i+1][j])*Uy1 + Uy[i+1][j]

                    Uy2 = 1.0
                    Ux2 = fabs(float(z1[i+1][j]))/(fabs(float(z1[i+1][j+1]))+fabs(float(z1[i+1][j])))
                    Uxx2 = (Ux[i+1][j+1] - Ux[i+1][j])*Ux2 + Ux[i+1][j]
                    Uyy2 = (Uy[i+1][j+1] - Uy[i+1][j])*Ux2 + Uy[i+1][j]

            if (mark_11[i][j+1] != 0.0):
                if (mark_12[i][j] != 0.0):
                    Ux1 = 1.0
                    Uy1 = fabs(float(z1[i][j+1]))/(fabs(float(z1[i][j+1]))+fabs(float(z1[i+1][j+1])))
                    Uxx1 = (Ux[i][j+1] - Ux[i+1][j+1])*Uy2 + Ux[i+1][j+1]
                    Uyy1 = (Uy[i][j+1] - Uy[i+1][j+1])*Uy2 + Uy[i+1][j+1]


                    Uy2 = 0.0
                    Ux2 = fabs(float(z1[i][j]))/(fabs(float(z1[i][j+1]))+fabs(float(z1[i][j])))
                    Uxx2 = (Ux[i][j+1] - Ux[i][j])*Ux2 + Ux[i][j]
                    Uyy2 = (Uy[i][j+1] - Uy[i][j])*Ux2 + Uy[i][j]

                if (mark_12[i+1][j] != 0.0):
                    Ux1 = 1.0
                    Uy1 = fabs(float(z1[i][j+1]))/(fabs(float(z1[i][j+1]))+fabs(float(z1[i+1][j+1])))
                    Uxx1 = (Ux[i][j+1] - Ux[i+1][j+1])*Uy2 + Ux[i+1][j+1]
                    Uyy1 = (Uy[i][j+1] - Uy[i+1][j+1])*Uy2 + Uy[i+1][j+1]

                    Uy2 = 1.0
                    Ux2 = fabs(float(z1[i+1][j]))/(fabs(float(z1[i+1][j+1]))+fabs(float(z1[i+1][j])))
                    Uxx2 = (Ux[i+1][j+1] - Ux[i+1][j])*Ux2 + Ux[i+1][j]
                    Uyy2 = (Uy[i+1][j+1] - Uy[i+1][j])*Ux2 + Uy[i+1][j]

            if (mark_12[i][j] != 0.0):
                if (mark_12[i+1][j] != 0.0):
                    Uy1 = 0.0
                    Ux1 = fabs(float(z1[i][j]))/(fabs(float(z1[i][j+1]))+fabs(float(z1[i][j])))
                    Uxx1 = (Ux[i][j+1] - Ux[i][j])*Ux2 + Ux[i][j]
                    Uyy1 = (Uy[i][j+1] - Uy[i][j])*Ux2 + Uy[i][j]

                    Uy2 = 1.0
                    Ux2 = fabs(float(z1[i+1][j]))/(fabs(float(z1[i+1][j+1]))+fabs(float(z1[i+1][j])))
                    Uxx2 = (Ux[i+1][j+1] - Ux[i+1][j])*Ux2 + Ux[i+1][j]
                    Uyy2 = (Uy[i+1][j+1] - Uy[i+1][j])*Ux2 + Uy[i+1][j]

            #------------------------------------------------------------------------

            if (mark_21[i][j] != 0.0):
                if (mark_21[i][j+1] != 0.0):
                    Qx1 = 0.0
                    Qy1 = fabs(float(z2[i][j]))/(fabs(float(z2[i][j]))+fabs(float(z2[i+1][j])))
                    Qxx1 = (Qx[i][j] - Qx[i+1][j])*Qy1 + Qx[i+1][j]
                    Qyy1 = (Qy[i][j] - Qy[i+1][j])*Qy1 + Qy[i+1][j]

                    Qx2 = 1.0
                    Qy2 = fabs(float(z2[i][j+1]))/(fabs(float(z2[i][j+1]))+fabs(float(z2[i+1][j+1])))
                    Qxx2 = (Qx[i][j+1] - Qx[i+1][j+1])*Qy2 + Qx[i+1][j+1]
                    Qyy2 = (Qy[i][j+1] - Qy[i+1][j+1])*Qy2 + Qy[i+1][j+1]

                if (mark_22[i][j] != 0.0):
                    Qx1 = 0.0
                    Qy1 = fabs(float(z2[i][j]))/(fabs(float(z2[i][j]))+fabs(float(z2[i+1][j])))
                    Qxx1 = (Qx[i][j] - Qx[i+1][j])*Qy1 + Qx[i+1][j]
                    Qyy1 = (Qy[i][j] - Qy[i+1][j])*Qy1 + Qy[i+1][j]

                    Qy2 = 0.0
                    Qx2 = fabs(float(z2[i][j]))/(fabs(float(z2[i][j+1]))+fabs(float(z2[i][j])))
                    Qxx2 = (Qx[i][j+1] - Qx[i][j])*Ux2 + Qx[i][j]
                    Qyy2 = (Qy[i][j+1] - Qy[i][j])*Ux2 + Qy[i][j]

                if (mark_22[i+1][j] != 0.0):
                    Qx1 = 0.0
                    Qy1 = fabs(float(z2[i][j]))/(fabs(float(z2[i][j]))+fabs(float(z2[i+1][j])))
                    Qxx1 = (Qx[i][j] - Qx[i+1][j])*Qy1 + Qx[i+1][j]
                    Qyy1 = (Qy[i][j] - Qy[i+1][j])*Qy1 + Qy[i+1][j]

                    Qy2 = 1.0
                    Qx2 = fabs(float(z2[i+1][j]))/(fabs(float(z2[i+1][j+1]))+fabs(float(z2[i+1][j])))
                    Qxx2 = (Qx[i+1][j+1] - Qx[i+1][j])*Qx2 + Qx[i+1][j]
                    Qyy2 = (Qy[i+1][j+1] - Qy[i+1][j])*Qx2 + Qy[i+1][j]

            if (mark_21[i][j+1] != 0.0):
                if (mark_22[i][j] != 0.0):
                    Qx1 = 1.0
                    Qy1 = fabs(float(z2[i][j+1]))/(fabs(float(z2[i][j+1]))+fabs(float(z2[i+1][j+1])))
                    Qxx1 = (Qx[i][j+1] - Qx[i+1][j+1])*Qy2 + Qx[i+1][j+1]
                    Qyy1 = (Qy[i][j+1] - Qy[i+1][j+1])*Qy2 + Qy[i+1][j+1]


                    Qy2 = 0.0
                    Qx2 = fabs(float(z2[i][j]))/(fabs(float(z2[i][j+1]))+fabs(float(z2[i][j])))
                    Qxx2 = (Qx[i][j+1] - Qx[i][j])*Ux2 + Qx[i][j]
                    Qyy2 = (Qy[i][j+1] - Qy[i][j])*Ux2 + Qy[i][j]

                if (mark_22[i+1][j] != 0.0):
                    Qx1 = 1.0
                    Qy1 = fabs(float(z2[i][j+1]))/(fabs(float(z2[i][j+1]))+fabs(float(z2[i+1][j+1])))
                    Qxx1 = (Qx[i][j+1] - Qx[i+1][j+1])*Qy2 + Qx[i+1][j+1]
                    Qyy1 = (Qy[i][j+1] - Qy[i+1][j+1])*Qy2 + Qy[i+1][j+1]

                    Qy2 = 1.0
                    Qx2 = fabs(float(z2[i+1][j]))/(fabs(float(z2[i+1][j+1]))+fabs(float(z2[i+1][j])))
                    Qxx2 = (Qx[i+1][j+1] - Qx[i+1][j])*Qx2 + Qx[i+1][j]
                    Qyy2 = (Qy[i+1][j+1] - Qy[i+1][j])*Qx2 + Qy[i+1][j]

            if (mark_22[i][j] != 0.0):
                if (mark_22[i+1][j] != 0.0):
                    Qy1 = 0.0
                    Qx1 = fabs(float(z2[i][j]))/(fabs(float(z2[i][j+1]))+fabs(float(z2[i][j])))
                    Qxx1 = (Qx[i][j+1] - Qx[i][j])*Ux2 + Qx[i][j]
                    Qyy1 = (Qy[i][j+1] - Qy[i][j])*Ux2 + Qy[i][j]

                    Qy2 = 1.0
                    Qx2 = fabs(float(z2[i+1][j]))/(fabs(float(z2[i+1][j+1]))+fabs(float(z2[i+1][j])))
                    Qxx2 = (Qx[i+1][j+1] - Qx[i+1][j])*Qx2 + Qx[i+1][j]
                    Qyy2 = (Qy[i+1][j+1] - Qy[i+1][j])*Qx2 + Qy[i+1][j]



            k1 = (Uy2 - Uy1)/(Ux2 - Ux1)
            k2 = (Qy2 - Qy1)/(Qx2 - Qx1)
        
            x_solve = (k1*(Ux1) - k2*(Qx1) + Qy1 - Uy1)/(k1 - k2)
            y_solve = k1*(x_solve - Ux1) + Uy1
            
            Uxx = (Uxx2 - Uxx1)/(Ux2 - Ux1)*x_solve + Uxx1
            Uyy = (Uyy2 - Uyy1)/(Uy2 - Uy1)*y_solve + Uyy1

            Qxx = (Qxx2 - Qxx1)/(Qx2 - Qx1)*x_solve + Qxx1
            Qyy = (Qyy2 - Qyy1)/(Qy2 - Qy1)*y_solve + Qyy1

            if (x_solve > -0.0) & (x_solve < 1.0):
                if (y_solve > -0.0) & (y_solve < 1.0):
                    answer_x.append(float(i)/N)
                    answer_y.append(float(j)/N)
                    answer_Ux.append(Uxx)
                    answer_Uy.append(Uyy)
                    answer_Qx.append(Qxx)
                    answer_Qy.append(Qyy)

#pcolor(x,y,sqrt(z1*z1 + z2*z2))

#plot(answer_x, answer_y, 'kx', ms=20)

#a1, b1 = genfromtxt('../data_mnk_f/lev1.dat').T
#plot(a1, b1, 'ko', ms=1)

#a2, b2 = genfromtxt('../data_mnk_f/lev2.dat').T
#plot(a2, b2, 'ko', ms=1)

#show()

file = open("pointsP.dat", 'w')
for k in range(0, size(answer_x)-1):
    file.write('{}    {}\n'.format(int(N*answer_x[k]), int(N*answer_y[k])))
file.close()

#-----------------------------------------

from numpy import roots, random
from math import atan, fabs


def cubic (Qx, Qy, Ux, Uy):
    
    a = Uy
    b = (Ux + 2*Qy)
    c = (2*Qx - Uy)
    d = -Ux
    
    det = -4*b*b*b*d + b*b*c*c -4*a*c*c*c + 18*a*b*c*d - 27*a*a*d*d
    
    if (det < 0):
        return 'c'
        
    
    if (det > 0):
        
        a = roots([a, b, c, d])
        a = a.real
        a = [atan(a[0]), atan(a[1]), atan(a[2])]
        
        U = [Ux*cos(a[0]) + Uy*sin(a[0]), Ux*cos(a[1]) + Uy*sin(a[1]), Ux*cos(a[2]) + Uy*sin(a[2])]
        rightU = [2*sin(a[0])*cos(a[0]), 2*sin(a[1])*cos(a[1]), 2*sin(a[2])*cos(a[2])]
    
        for i in range(0, 3):
            if (U[i] * rightU[i] < 0):
                a[i] = a[i] + pi
        
        a = sorted(a)
        a = [a[0] - a[0], a[1] - a[0], a[2] - a[0]]
        
        #print a
        
        if (a[2] > pi):
            return 'a'
        else:
            return 'b'


#-----------------------------------------

type = []

number = [0, 0, 0]

for k in range(0, size(answer_x)-1):
    type.append(cubic(answer_Qx[k], answer_Qy[k], answer_Ux[k], answer_Uy[k]))
    if (cubic(answer_Qx[k], answer_Qy[k], answer_Ux[k], answer_Uy[k]) == 'a'):
        number[0] = number[0] + 1
    if (cubic(answer_Qx[k], answer_Qy[k], answer_Ux[k], answer_Uy[k]) == 'b'):
        number[1] = number[1] + 1
    if (cubic(answer_Qx[k], answer_Qy[k], answer_Ux[k], answer_Uy[k]) == 'c'):
        number[2] = number[2] + 1

#print np.array(number) /float(number[0] + number[1] + number[2])
file = open("type_A.dat", 'a')
file.write('{}\n'.format(number[0]/float(number[0] + number[1] + number[2])))
file = open("type_B.dat", 'a')
file.write('{}\n'.format(number[1]/float(number[0] + number[1] + number[2])))
file = open("type_C.dat", 'a')
file.write('{}\n'.format(number[2]/float(number[0] + number[1] + number[2])))


file = open("points_typesP.dat", 'w')
for k in range(0, size(answer_x)-1):
    file.write('{}    {}    {}\n'.format(int(N*answer_x[k]), int(N*answer_y[k]), type[k]))

#------------------------------------------

