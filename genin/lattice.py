import sys
import numpy
import os.path
import logging
kennyloggins = logging.getLogger('lattice')

class Lattice:
    @staticmethod
    def fromfile(fname, m):
        kennyloggins.debug("Loading lattice from file " + os.path.abspath(fname))
        latfile = open(fname,'r')
        f = latfile.readlines()
        latfile.close()
        #Some lattice files will have the node they ran on printed at the top; bypass 1st line by trying again
        try: n = int(f.pop(0))
        except ValueError:
            n = int(f.pop(0))

        h = numpy.zeros((n,n))
        for pl,item in enumerate(f):
            s = item.split()
            if len(s) == 3:
                i,j = int(s[0]), int(s[1])
                v = float(s[2])
                h[i,j] = v
            if len(s) != 3:
                f = f[pl:]
                break

        V_check = f[-1].split()
        if V_check[0] == 'U':
            V = float(V_check[1])
        elif len(V_check) == 5:
            V = numpy.zeros((n,n,n,n))
            for item in f:
                s = item.split()
                i,j,k,l = [int(ind) for ind in s[:-1]]
                v = float(s[4])
                V[i,j,k,l] = v
        return Lattice(n,h,V,m)

    def __init__(s,n,h,V,m):
        s.m = m
        s.n = n
        s.h = h
        s.V = V
        if not (s.n == len(s.h)):
            kennyloggins.error("Number of sites and integrals misaligned")
            import sys; sys.exit(1)

    def filling(s):
        return float(s.m) / s.n


