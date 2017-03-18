#Code doesn't work

import random
from math import *

points = float(254)
modes = float(3)
seed = 5

file = open('foo.dat', 'w')
random.seed()

def f(k1, k2, x, y):
    return cos(k1*x+k2*y) + sin(k1*x+k2*y)


for j1 in range(0, int(points)):
    for j2 in range (0, int(points)):
        fx_ = 0.0
        x = 2*pi*j1/points
        y = 2*pi*j2/points
        fx_ = random.random()*f(1, 1, x, y) + random.random()*f(1, 2, x, y) + random.random()*f(1, 3, x, y) + random.random()*f(2, 2, x, y) + random.random()*f(2, 3, x, y) + random.random()*f(3, 3, x, y)
            
        file.write(repr(fx_) + '  ' + repr(j1/points) + '   ' + repr(j2/points) + '\n')
