import sys
import numpy
sys.path.append('..')

import solvers
import config
import logging
logging.basicConfig()

n = 10
h = numpy.diag(numpy.ones(n-1),1)
h[0,-1] = 1
h+=h.T
h*=-1

V = numpy.zeros((n,n,n,n))
for i in range(n):
    V[i,i,i,i] = 8

fci = solvers.TroyFCI()
fci.solve(h,V)
print fci.get_energy()

dh = numpy.zeros((n,n))
dh[0,0] = 1e-6

config.solvers['recycle_guess'] = True
fci.solve(h+dh,V)
print fci.get_energy()
PA = fci.get_P()
print PA

config.solvers['recycle_guess'] = False
fci.solve(h+dh,V)
print fci.get_energy()
PB = fci.get_P()
print PB

print PA-PB
print numpy.linalg.norm(PA-PB)
