import numpy
import scipy.io
import subprocess
import uuid
import os
import shutil

n = 5
ns = 1
#h = numpy.random.randn(2*n,2*n)
#V = numpy.random.randn(2*n,2*n,2*n,2*n)
h = numpy.diag(numpy.ones(2*n-1),1)
h[0,-1] = 1
h += h.T
h *= -1

V = numpy.zeros((2*n, 2*n, 2*n, 2*n))
for i in range(2*n):
    V[i,i,i,i] = 8

infname = os.path.join('/dev/shm/', str(uuid.uuid4()))
outfname = os.path.join('/dev/shm/',str(uuid.uuid4()))
f = scipy.io.FortranFile(infname,'w')
f.write_record(h.T)
f.write_record(V.T)
f.close()
#p = subprocess.Popen(['./fci', str(n), str(ns)], stdout=subprocess.PIPE, stdin=subprocess.PIPE)

import sh
sh.Command('./fci')(n,ns,infname,outfname)

f = scipy.io.FortranFile(outfname,'r')
Pout = f.read_reals().reshape(ns,ns,2*n,2*n)
OnTopout = f.read_reals().reshape(ns,ns,2*n)
Eout = f.read_reals()
f.close()

Pout = Pout.T
OnTopout = OnTopout.T

sh.rm(infname,outfname)

print Eout
print Pout[:,:,0,0]

