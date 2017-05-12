### Test transformation of 2 e- integrals
import sys
sys.path.append('..')

import numpy
import logging
logging.basicConfig()

import xform
import testutil

U = 4
n = 10
ni = 3
Q,R = numpy.linalg.qr(numpy.random.randn(n,n))
T = Q[:,:ni]
Vfull = numpy.zeros((n,n,n,n))
Vline = numpy.zeros(n)
for i in range(n):
    Vfull[i,i,i,i] = U
    Vline[i] = U

V1 = xform.two(T,Vfull)
V2 = xform.two(T,Vline)
V3 = xform.twoU(T,U)

testutil.print_result( (V1 == V2).all() and (V2 == V3).all() )
