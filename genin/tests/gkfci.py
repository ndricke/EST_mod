import sys
import numpy
sys.path.append('..')

import solvers
import hf
import logging
import lattice
logging.basicConfig()


l = lattice.Lattice(10,1,5,8)
h,U = l.get_operators()

V = numpy.zeros((10,10,10,10))
for i in range(10):
    V[i,i,i,i] = U
gkfci = solvers.GKFCI()
gkfci.solve(h,V,10)
print gkfci.energy

