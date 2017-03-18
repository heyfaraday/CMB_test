#!/usr/bin/python
# -*- coding: UTF-8 -*-
 
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import argparse
import subprocess
 
def createParser ():
    parser = argparse.ArgumentParser()
    parser.add_argument ('-c', '--count', default=1, type=int)
 
    return parser
 
 

if (__name__ == '__main__'):
    parser = createParser()
    namespace = parser.parse_args(sys.argv[1:])
 
    print (namespace)
 
for _ in range (namespace.count):
    cmd = 'python ifft.py'
    p = subprocess.Popen(cmd, shell = True)
    p.wait()
    
    cmd = './a.out'
    p = subprocess.Popen(cmd, shell = True)
    p.wait()
        
    cmd = 'python searchP.py'
    p = subprocess.Popen(cmd, shell = True)
    p.wait()
    
    cmd = 'python mean.py'
    p = subprocess.Popen(cmd, shell = True)
    p.wait()
    
    cmd = 'python mapP.py'
    p = subprocess.Popen(cmd, shell = True)
    p.wait()