import sys
sys.path.append('..')

import numpy
import logging
logging.basicConfig(level=logging.DEBUG)

import al

def f(x):
    return x[0] + x[1], numpy.array([x[0]**2 + x[1]**2 - 2 ])

x = al.al(f,numpy.array([0,0]))
print x
print f(x)

