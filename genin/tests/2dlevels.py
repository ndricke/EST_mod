import numpy
import argparse
import matplotlib.pyplot as plt
import logging
logging.basicConfig()

import sys
sys.path.append('..')
import lattice

parser = argparse.ArgumentParser()
parser.add_argument('-nx', help='Number of sites', type=int)
parser.add_argument('-ny', help='Number of sites', type=int)
parser.add_argument('--pbcx', help='Periodic boundary conditions', type=int, default=1)
parser.add_argument('--pbcy', help='Periodic boundary conditions', type=int, default=-1)

args = parser.parse_args()

l = lattice.Lattice(args.nx, args.ny, 0, 0, args.pbcx, args.pbcy)

h,U = l.get_operators()

w,v = numpy.linalg.eigh(h)

plt.hist(w, bins=args.nx*args.ny*10)
plt.show()
for m in range(7,args.nx*args.ny/2+1):
    if abs(w[m-1] - w[m]) > 1e-2:
        print m
for i in w: print i
