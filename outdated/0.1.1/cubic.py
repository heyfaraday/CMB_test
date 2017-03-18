from numpy import roots
from math import atan, fabs

def cubic (Qx, Qy, Ux, Uy):
    
    a = Uy
    b = (Ux + 2*Qy)
    c = (2*Qx - Uy)
    d = -Ux
    
    det = -4*b*b*b*d + b*b*c*c -4*a*c*c*c + 18*a*b*c*d - 27*a*a*d*d
    
    print det
    
    if (det < 0):
        answer = 'c'
        #print answer
        
    if (det > 0):
        a = roots([a, b, c, d])
        a = a.real
        #print a
        a = atan(a[0])/(fabs(atan(a[0]))) + atan(a[1])/(fabs(atan(a[1]))) + atan(a[2])/(fabs(atan(a[2])))
        
        if (a == 3.0):
            answer = 'b'
        if (a == 1.0):
            answer = 'b'
        if (a == -1.0):
            answer = 'a'
        #print answer
        
    return answer