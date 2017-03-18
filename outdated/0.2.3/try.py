import numpy as np
from math import sin, cos, tan, pi, sqrt, factorial, fabs, acos

from numpy import complex128, float64
#import time
import pyfftw
from pyfftw.pyfftw import FFTW

Lmax = 6

a_coef = np.random.normal(size = (Lmax+1, Lmax+1))
b_coef = np.random.normal(size = (Lmax+1, Lmax+1))

#a_coef = np.ones((Lmax+1, Lmax+1))
#b_coef = np.ones((Lmax+1, Lmax+1))

a_coef[0][0] = 0.0

for m in xrange(0, Lmax+1):
    for l in xrange(0, m):
        a_coef[m][l] = 0.0     
        
for m in xrange(0, Lmax+1):
    for l in xrange(0, m):
        b_coef[m][l] = 0.0
        
for l in xrange(0, Lmax+1):
        b_coef[0][l] = 0.0

N = 512

# P_

P_ = np.zeros((N/2+1, Lmax+4, Lmax+4))

for j in xrange(1, N/2):
    
    teta = 2*pi*j/float(N)
    
    P_[j][0][0] = 1/sqrt(4*pi)

    for m in xrange(1, Lmax+1):
            P_[j][m][m] = P_[j][m-1][m-1]*(-sin(teta))*sqrt(2*m+3)/sqrt(2*m+2)
    
    for m in xrange(0, Lmax):
            P_[j][m][m+1] = P_[j][m][m]*cos(teta)*sqrt(2*m+3)
    
    for m in xrange(0, Lmax-1):
            for l in xrange(m+2, Lmax+1):
                P_[j][m][l] = sqrt((2*l+1)*(l-1-m))/sqrt(l**2-m**2)*(cos(teta)*sqrt(2*l-1)/sqrt(l-1-m)*P_[j][m][l-1] - sqrt(l+m-1)/sqrt(2*l-3)*P_[j][m][l-2])

x = np.zeros((N, N/2))
y = np.zeros((N, N/2))

for j in xrange(1,N/2):
	for i in xrange(0, N):
        phi = pi*i*2/float(N)
           
        field[i][j] = T[i]
            
        x[i][j] = (i-N/2)*2/float(N)*pi
        y[i][j] = 2*pi*j/float(N)- pi/2*(N/4)*4/float(N)

field = np.zeros((N, N/2))
                
for j in xrange(1, N/2):
    
    teta = 2*pi*j/float(N)
                
    F = complex128(np.zeros((N+1)))
    F_ = complex128(np.zeros((N+1)))   
        
    func1 = 0.0
    func2 = 0.0
        
    for m in xrange(0, Lmax+1):
        for l in xrange(m, Lmax+1):
            func1 = func1 + a_coef[m][l]*P_[j][m][l]
            func2 = func2 + b_coef[m][l]*P_[j][m][l]
    
        F[m] = func1
        F_[m] = func2
                
        func1 = 0.0
        func2 = 0.0
            
    T = np.real(pyfftw.interfaces.numpy_fft.fft(F)) + np.imag(pyfftw.interfaces.numpy_fft.fft(F_))
        
    for i in xrange(0, N):
        phi = pi*i*2/float(N)
           
        field[i][j] = T[i]
            
            

