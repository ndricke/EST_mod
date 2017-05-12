import numpy
import itertools
import logging

if __name__=="__main__": logging.basicConfig()
kennyloggins = logging.getLogger('hf')

class HF:
    def __init__(s,h,V,m,diis=True,tol=1e-10,rho_guess=None):
        s.h = h
        s.V = V
        s.m = m
        s.diis = diis
        s.tol = tol
        s.count = 0

        n = s.h.shape[0] #Num rows is the dimension of the density matrix to be
        if type(rho_guess) == numpy.ndarray:
            s.rho = rho_guess
            if s.rho.shape[0] != n: raise ValueError("Input HF guess has wrong shape")
        else:
            s.rho = numpy.eye(n)

        s.C = None
        s.e = None

    def _do_hf(s):
        n = s.h.shape[0]
        ## Assume Hubbard interaction
        if type(s.V) == float:
            F = s.h + numpy.diag(numpy.diag(s.rho)*s.V)
        elif len(s.V.shape)==4:
            F = s.h.copy()
            ri = range(len(s.h))
            for mu,nu,la,si in itertools.product(ri,ri,ri,ri):
                F[mu,nu] += (2*s.V[mu,si,nu,la] - s.V[mu,si,la,nu])*s.rho[la,si]
        else:
            print "Unknown V shape"
            import sys
            sys.exit(1)

        converged = False
        niter = 0
        errs = []
        Fs = []
        while not converged:
            s.count += 1; print s.count
            w,C = numpy.linalg.eigh(F)
            
            # Check for HOMO-LUMO degeneracy
            homo = w[s.m-1]
            lumo = w[s.m]
            if abs(lumo-homo) < 1e-4:
                kennyloggins.warning("DANGER ZONE: HOMO and LUMO are degenerate")

            Cocc = C[:,:s.m]
            rho_new = Cocc.dot(Cocc.T)

            ## Assume Hubbard interaction
            if type(s.V) == float:
                F_new = s.h + numpy.diag(numpy.diag(rho_new)*s.V)
            elif len(s.V.shape)==4:
                F_new = s.h.copy()
                ri = range(len(s.h))
                for mu,nu,la,si in itertools.product(ri,ri,ri,ri):
                    F_new[mu,nu] += (2*s.V[mu,si,nu,la] - s.V[mu,si,la,nu])*rho_new[la,si]
            else:
                print "Unknown V shape"
                import sys
                sys.exit(1)
            F_new_MO = C.T.dot(F_new).dot(C)
            #err = F_new_MO[:m,m:]
            err = F_new.dot(rho_new) - rho_new.dot(F_new)


            errs.append(err)
            Fs.append(F_new)
            if s.diis:
                F = HF._next_diis(errs, Fs)
            else:
                s.rho = rho_new
                F = F_new

            converged = numpy.linalg.norm(err)<s.tol*n
            niter += 1
            #if niter > 3:
            #    asdflkj
            emat = Cocc.T.dot(s.h+F).dot(Cocc)
            kennyloggins.debug("Levels: %s" % str([emat[i,i] for i in range(len(emat))]))
        emat = Cocc.T.dot(s.h+F).dot(Cocc)
        e = 0
        for i in range(s.m):
            e += emat[i,i]
        s.C = C
        s.e = e

    @staticmethod
    def _next_diis(errs, Fs):
        n = len(errs)
        B = numpy.zeros((n,n))
        for i,j in itertools.product(range(n), range(n)):
            B[i,j] = numpy.dot(errs[i].ravel(), errs[j].ravel())
        A = numpy.zeros((n+1, n+1))
        A[:n,:n] = B
        A[n,:] = -1
        A[:,n] = -1
        A[n,n] = 0
        b = numpy.zeros(n+1)
        b[n] = -1
        try:
            x = numpy.linalg.solve(A,b)
        except (numpy.linalg.linalg.LinAlgError):
            print "lin solver fails! Using pinv..."
            P = numpy.linalg.pinv(A)
            x = P.dot(b)
        w = x[:n]

        F = numpy.zeros(Fs[0].shape)
        for i in range(n):
            F += w[i] * Fs[i]
        return F

    def get_Cocc(s):
        if s.C is None: s._do_hf()
        return s.C[:,:s.m]

    def get_C(s):
        if s.C is None: s._do_hf()
        return s.C

    def get_energy(s):
        if s.e is None: s._do_hf()
        return s.e

    def get_rho(s):
        C = s.get_Cocc()
        return C.dot(C.T)

    def get_P2(s):
        P1 = s.get_rho()
        s.P2 = s.pdm2(P1)
        return s.P2

    def pdm2(s,P1):
        n = P1.shape[0]
        P2 = numpy.zeros((n,n,n,n))
        sz = range(n)
        for i,j,k,l in itertools.product(sz,sz,sz,sz):
            P2[i,j,k,l] = P1[k,i]*P1[j,l]
        return P2

if __name__ == "__main__":
    import lattice
    import sys
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', help='Number of electron pairs', type=int, required=True)
    parser.add_argument('-l', help='Lattice file', type=str, required=True)
    parser.add_argument('-o', help='Output file', type=str, required=True)
    parser.add_argument('-guess', help='HF Guess', type=str)
    args = parser.parse_args()

    rho = numpy.genfromtxt(args.guess)

    lat = lattice.Lattice.fromfile(args.l,args.m)
    if args.guess != None: fock = HF(lat.h,lat.V,lat.m,rho_guess = rho) 
    else: fock = HF(lat.h,lat.V,lat.m)
    rho = fock.get_rho()
    numpy.savetxt(args.o,rho)


#    numpy.set_printoptions(precision=3, suppress = True)
#    n = 480
#    m = 311
#    U = 8.0
#    h = numpy.diag(numpy.ones(n-1),1)
#    h[0,-1] = 1
#    h += h.T
#    h*=-1
#    for i in range(0,n,2):
#        h[i,i] = 1
#
#    hf = HF(h,U,m)
#    psi = hf.get_Cocc()
#    print psi.dot(psi.T)*2
#    print hf.get_energy()
