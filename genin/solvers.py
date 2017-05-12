import numpy
import os
import sh
import tempfile
import shutil
import logging
if __name__ == "__main__": logging.basicConfig()
kennyloggins = logging.getLogger('solver')
from troyfci import TroyFCI

# Deprecated
def takashi_fci(h,V):
    n = h.shape[0]
    f = open("FCIINP_HUB2", 'w')
    print >> f, " &FCI NORB %d NELEC %d MS2 0" % (n,n)
    print >> f, "  ORBSYM=",
    for i in range(n):
        print >> f, "1,"
    print >>f
    print >> f, "  ISYM=1"
    print >> f, " &END"

    for i in range(n):
        for j in range(n):
            print >> f, "%f %d %d 0 0" % (h[i,j], i+1, j+1)

    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    print >> f, "%f %d %d %d %d" % (V[i,j,k,l], i+1, j+1, k+1, l+1)
    
    f.close()

    os.system('/home/welborn/work/dmet/newshiny/matt/takashi_fci/fci')

    f = open('1RDM')
    f.readline()
    p_flat = numpy.array(map(float, f.readlines()))
    p = numpy.reshape(p_flat, (n,n))
    f.close()

    return p

class GKFCI: # Subclass?
    def __init__(s,usetmp=True,binary="source"):
        s.usetmp = usetmp

        if binary == "source": # load binary from source directory
            binary = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gkfci')
        s.fci = sh.Command(binary)

    def solve(s,h,V,nelec):
        if s.usetmp:
            tdir = tempfile.mkdtemp(dir="/scratch/welborn")
            pwd = os.getcwd()
            os.chdir(tdir)

        n = h.shape[0]
        f = open("FCIINP_HUB2", 'w')
        print >> f, " &FCI NORB= %d, NELEC= %d, MS2= 0," % (n,nelec)
        print >> f, "  ORBSYM="+"1,"*n
        print >>f
        print >> f, "  ISYM=1"
        print >> f, " &END"


        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        # change from physicists' notation to chemists'
                        print >> f, "%30.20f %d %d %d %d" % (V[i,k,j,l], i+1, j+1, k+1, l+1)

        for i in range(n):
            for j in range(n):
                print >> f, "%30.20f %d %d 0 0" % (h[i,j], i+1, j+1)
        
        f.close()

        try:
            out = s.fci('-v', "1e-11",'--basis', 'Input','--save-rdm1', '1rdm', '--save-rdm2', '2rdm', '--save-rdm2s', '2rdms', '--ptrace', n/2 , '--diis-block-size', '20000', 'FCIINP_HUB2')
        except sh.ErrorReturnCode_1:
            kennyloggins.error("GKFCI failed. Copying input to working directory...")
            shutil.copy('FCIINP_HUB2', os.path.join(pwd,"FCIINP.err"))
            import sys
            sys.exit(1)
        s.energy = numpy.nan
        s.ptrace = numpy.nan
        f = open(os.path.join(pwd,'fci.log'), 'w')
        for line in out:
            print >> f, line,
            if "Time for main loop:" in line:
                kennyloggins.debug("GKFCI time: "+ line.split()[-2])
            if "FCI STATE 1 pTraceSys" in line:
                s.ptrace = float(line.split()[-1])
            if "FCI STATE 1 ENERGY" in line:
                s.energy = float(line.split()[-1])
                break
        else:
            kennyloggins.warning("Failed to get energy from FCI")
        f.close()

        f = open('1rdm')
        line = f.readline()
        assert(int(line.split()[-1]) == n)
        s.P = numpy.array([map(float,line.split()) for line in f])
        f.close()

        f = open('2rdm')
        line = f.readline()
        s.P2 = numpy.array([map(float,line.split()) for line in f]).reshape(n,n,n,n)
        f.close()

        f = open('2rdms')
        line = f.readline()
        s.P2s = numpy.array([map(float,line.split()) for line in f]).reshape(n,n,n,n)
        f.close()

        if s.usetmp:
            os.chdir(pwd)
            shutil.rmtree(tdir)
    def get_energy(s, state=0):
        if state != 0:
            kennyloggins.error("GKFCI only does ground states")
            import sys
            sys.exit(-1)
        return s.energy
    def get_energies(s):
        return [s.energy]
    def get_P(s, state1=0, state2=0):
        if state1 != 0 and state2 != 0:
            kennyloggins.error("GKFCI only does ground states")
            import sys
            sys.exit(-1)
        return s.P
    def get_Ps(s):
        return [[s.P]]
    def get_OnTop(s, state1=0, state2=0):
        if state1 != 0 and state2 != 0:
            kennyloggins.error("GKFCI only does ground states")
            import sys
            sys.exit(-1)
        n = s.P.shape[0]
        OnTop = numpy.zeros(n)
        for i in range(n):
            OnTop[i] = s.P2[i,i,i,i] - s.P2s[i,i,i,i]
        return OnTop
    def get_OnTops(s):
        raise NotImplementedError


if __name__ == "__main__":
    # Test Troy's and GK's fci on the hubbard model

    n = 8
    h = numpy.diag(numpy.ones(n-1),1)
    h[0,-1] = 1
    h += h.T
    h *= -1

    U = numpy.zeros((n,n,n,n))
    for i in range(n):
        U[i,i,i,i] = 8

    gkfci = GKFCI()
    gkfci.solve(h,U,n)
    print "GKFCI"
    print gkfci.get_energy()
    print gkfci.get_P()
    print gkfci.get_OnTop()

    print

    print "Troy FCI"
    tfci = TroyFCI()
    tfci.solve(h,U,n)
    print tfci.get_energy()
    print tfci.get_P()
    print tfci.get_OnTop()




