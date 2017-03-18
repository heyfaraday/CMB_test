# Code works stable with environment.yml 

from numpy import *
from pylab import *
import pandas as pd

# a,b = genfromtxt('plot.dat').T

N = 254 
M = 254

xlim([0,1])
ylim([0,1])
 
# axes().set_aspect('equal', 'datalim')


x,y,z = genfromtxt('../data_map_f/f1.dat').T
x = x.reshape(N,M)
y = y.reshape(N,M)
z = z.reshape(N,M)

x0,y0,z0 = genfromtxt('../data_map_f/f2.dat').T
x0 = x.reshape(N,M)
y0 = y.reshape(N,M)
z0 = z.reshape(N,M)

pcolor(x, y, np.sqrt(z0*z0 + z*z))

# We are search the zero points on both zero contour lines
a, b = genfromtxt('../data_mnk_f/lev1.dat').T 
a0, b0 = genfromtxt('../data_mnk_f/lev2.dat').T
plot(a, b,'ko',ms=1)
plot(a0, b0,'ko',ms=1)

list1 = []

for i in range(0, 20000):
    for j in range(0, 20000):
        if ((a[i]-a0[j])**2 + (b[i] - b0[j])**2 < 0.000000005) & ((a[i]-a0[j])**2 + (b[i] - b0[j])**2 > 0.0000000005):
            print i,' ', j
            list1.append(i)

plot(a[list1], b[list1],'kx',ms=20)

#xmax,ymax = genfromtxt('max.dat').T
#plot(xmax,ymax,'ko',ms=1)

#xmin,ymin = genfromtxt('min.dat').T
#plot(xmin,ymin,'ko',ms=4)

#xsad,ysad = genfromtxt('sad.dat').T
#plot(xsad,ysad,'kx',ms=4)


#print 'Enter filename or nothing to show:'
#f = raw_input()
#if raw_input() != '':
#    savefig(f)
#else:

show()
#savefig('fig.png', dpi=1500)
