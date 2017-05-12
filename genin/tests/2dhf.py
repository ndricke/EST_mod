import sys
sys.path.append('..')

import lattice
import hf
import logging
logging.basicConfig()

nx = 16
ny = 16
m = nx*ny/2

print "Trying %d by %d w/ %d electrons..." % (nx,ny,m)
l = lattice.Lattice(nx,ny,m,8)
h,U = l.get_operators()
fock = hf.HF(h,U,m)
print fock.get_energy()
print "Done"
print


nx = 16
ny = 16
m = 52

print "Trying %d by %d w/ %d electrons..." % (nx,ny,m)
l = lattice.Lattice(nx,ny,m,8)
h,U = l.get_operators()
fock = hf.HF(h,U,m)
print fock.get_energy()
print "Done"
print
