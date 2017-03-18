from numpy import *
from pylab import *

A = genfromtxt('type_A.dat').T
B = genfromtxt('type_B.dat').T
C = genfromtxt('type_C.dat').T
print np.mean(A),  np.mean(B), np.mean(C)