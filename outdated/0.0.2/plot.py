# Code works stable with environment.yml 

from numpy import *
from pylab import *
import pandas as pd

N = 254
M = 254

xlim([0,1])
ylim([0,1])

x1,y1,z1 = genfromtxt('foo.dat').T
x2,y2,z2 = genfromtxt('foo_sort.dat').T

x1 = x1.reshape(N,M)
y1 = y1.reshape(N,M)
z1 = z1.reshape(N,M)

plot(x2[0:100],y2[0:100],'kx',ms=20)

pcolor(x1, y1, z1)

show()

