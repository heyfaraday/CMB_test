# Code works stable with environment.yml 

from numpy import *
from pylab import *
import pandas as pd

x1,y1,z1 = genfromtxt('../data_map_f/f1.dat').T
x2,y2,z2 = genfromtxt('../data_map_f/f2.dat').T

z = sqrt(z1*z1 + z2*z2)

d = np.array([x1, y1, z])
df = pd.DataFrame({'x': x1, 'y': y1, 'z': z})

df = df.sort_values(by=['z'], axis=0, ascending=True, inplace=False, kind='quicksort', na_position='last')

np.savetxt("foo.dat", d.T, delimiter="    ")
np.savetxt("foo_sort.dat", df, delimiter="    ")
