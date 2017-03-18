#work with new2 environment
from numpy import *
from pylab import *
import pandas as pd
#a,b = genfromtxt('plot.dat').T

N=254
M=254

xlim([0,1])
ylim([0,1])

#axes().set_aspect('equal', 'datalim')


x,y,z = genfromtxt('../data/f1.dat').T
x = x.reshape(N,M)
y = y.reshape(N,M)
z = z.reshape(N,M)

x0,y0,z0 = genfromtxt('../data/f1.dat').T
x0 = x.reshape(N,M)
y0 = y.reshape(N,M)
z0 = z.reshape(N,M)

pcolor(x, y,np.sqrt(z0*z0 + z*z))

a,b = genfromtxt('../data/lev1.dat').T
a0,b0 = genfromtxt('../data/lev2.dat').T
plot(a, b,'ko',ms=1)
plot(a0, b0,'ko',ms=1)


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