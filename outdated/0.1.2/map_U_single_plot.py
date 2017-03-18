from numpy import *
from pylab import *
#a,b = genfromtxt('plot.dat').T

xlim([0,1])
ylim([0,1])

#axes().set_aspect('equal', 'datalim')

x,y,z = genfromtxt('f2.dat').T
N = sqrt(x.shape[0])
x = x.reshape(N,N)
y = y.reshape(N,N)
z = z.reshape(N,N)
pcolor(x,y,z)

#a,b = genfromtxt('../data/lev.dat').T
#plot(a,b,'ko',ms=1)

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
