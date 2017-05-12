
import numpy
import scipy.io
import scipy.misc
import tempfile
import os
import shutil
import sh
import logging
if __name__ == "__main__": logging.basicConfig()
kennyloggins = logging.getLogger('solver')

import config

class TroyFCI: # Subclass?
    def __init__(s,usetmp=True,binary="source"):
        s.usetmp = usetmp

        if binary == "source": # load binary from source directory
            binary = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'troyci/fci')
        s.fci = sh.Command(binary)

        hred_binary = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'troytest/hred')
        s.hred = sh.Command(hred_binary)

    def solve(s,h,V,hnat=None,Vnat=None,nelec=None, ns=1):
        if s.usetmp:
            tdir = tempfile.mkdtemp(dir="/dev/shm")
            pwd = os.getcwd()
            os.chdir(tdir)
        else: tdir = "./"

        n = h.shape[0]

        if nelec is None: nelec = n//2
        if Vnat == None: Vnat = V
        if hnat == None: hnat = h

        infname = os.path.join(tdir, 'in.bin')
        outfname = os.path.join(tdir, 'out.bin')
        f = scipy.io.FortranFile(infname,'w')
        f.write_record(h.T)
        f.write_record(V.T)
        f.write_record(hnat.T)
        f.write_record(Vnat.T)
        if config.solvers['recycle_guess']:
            f.write_record(s.psi)
        f.close()

        s.fci(False,n,nelec,ns,config.solvers['recycle_guess'])

        f = scipy.io.FortranFile(outfname, 'r')
        Pout = f.read_reals().reshape(ns,ns,n,n).T
        P2out = f.read_reals().reshape(n,n,n,n).T
        Eout = f.read_reals()
        Psiout = f.read_reals() # Perhaps we don't really need to modify anything about psi?
        H1s = f.read_reals().reshape(4,4).T
        
        f.close()


        if s.usetmp:
            os.chdir(pwd)
            shutil.rmtree(tdir)

        s.P = Pout * 2.
        s.P2 = P2out * 4.
        s.energy = Eout
        s.psi = Psiout #+ numpy.random.randn(len(Psiout))*.0001
        s.H_red = H1s

    def calcHred(s,m,h,V,X1,X2,X3,X4,):
        n = h.shape[0]
        if s.usetmp:
            tdir = tempfile.mkdtemp(dir="/dev/shm")
            pwd = os.getcwd()
            os.chdir(tdir)
        else: tdir = "./"
        infname = os.path.join(tdir, 'in.bin')
        outfname = os.path.join(tdir, 'out.bin')
        f = scipy.io.FortranFile(infname,'w')
        f.write_record(h.T)
        f.write_record(V.T)
        f.write_record(X1.T)
        f.write_record(X2.T)
        f.write_record(X3.T)
        f.write_record(X4.T)
        f.close()

        s.hred(n,m)

        f = scipy.io.FortranFile(outfname, 'r')
        H1s = f.read_reals().reshape(4,4).T
#        X1 = f.read_reals().reshape(2,2).T
#        X2 = f.read_reals().reshape(2,2).T
#        X1 = f.read_reals()
#        X2 = f.read_reals().reshape(2,2,2,2).T
        f.close()

#        print X1
#        print X2

        if s.usetmp:
            os.chdir(pwd)
            shutil.rmtree(tdir)

        return H1s

    def calcPG(s,n,m,psi):
        if s.usetmp:
            tdir = tempfile.mkdtemp(dir="/dev/shm")
            pwd = os.getcwd()
            os.chdir(tdir)
        else: tdir = "./"
        infname = os.path.join(tdir, 'in.bin')
        outfname = os.path.join(tdir, 'out.bin')
        f = scipy.io.FortranFile(infname,'w')
        f.write_record(psi.T)
        f.close()

        s.fci(True,n,m,1)

        f = scipy.io.FortranFile(outfname, 'r')
        Pout = f.read_reals().reshape(n,n).T
        P2out = f.read_reals().reshape(n,n,n,n).T
        X = f.read_reals()
        f.close()

        if s.usetmp:
            os.chdir(pwd)
            shutil.rmtree(tdir)

        return Pout*2., P2out*4.

    def get_energy(s, state=0):
        return s.energy[state]
    def get_energies(s):
        return s.energy
    def get_P(s, state1=0, state2=0):
        return s.P[:,:,state1,state2]
    def get_Ps(s):
        return s.P
    def get_P2(s):
        return s.P2


if __name__ == "__main__":
    import lattice
    import sys

    lat_name = sys.argv[1]
    m = sys.argv[2]
    ns = 1

    l = lattice.Lattice.fromfile(lat_name,m)

    fci = TroyFCI(usetmp=False)
    fci.solve(l.h,l.V,nelec=l.m,ns=ns)

#    numpy.set_printoptions(precision=3, suppress=True)

    n = l.h.shape[0]
    print fci.energy/n












