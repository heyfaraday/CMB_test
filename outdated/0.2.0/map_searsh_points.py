from numpy import *
from pylab import *
#a,b = genfromtxt('plot.dat').T

N=254
M=254

xlim([0,1])
ylim([0,1])

#axes().set_aspect('equal', 'datalim')


z,x,y = genfromtxt('P.dat').T

x = x.reshape(N,M)
y = y.reshape(N,M)
z = z.reshape(N,M)

pcolormesh(x,y,z)

#a,b = genfromtxt('levP.dat').T
#plot(a,b,'ko',ms=1)

#xmax,ymax,fmax = genfromtxt('maxP.dat').T
#plot(xmax,ymax,'ko',ms=1)

xmin,ymin,fmin = genfromtxt('minP.dat').T
plot(xmin,ymin,'ko',ms=4)

xmin,ymin = genfromtxt('pointsP.dat').T
plot(xmin/float(N), ymin/float(N),'kx',ms=10)

#xsad,ysad = genfromtxt('sadP.dat').T
#plot(xsad,ysad,'kx',ms=4)


#print 'Enter filename or nothing to show:'
#f = raw_input()
#if raw_input() != '':
#    savefig(f)
#else:

show()

#savefig('fig.png', dpi=1500)
