from pylab import *
from numpy import *

x,y,z1 = genfromtxt('maxU.dat').T
plt.hist(z1)

x,y,z2 = genfromtxt('minU.dat').T
plt.hist(z2)
plt.show()
