#!/usr/bin/env python
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
        if rho_guess == None: s.rho = numpy.eye(n)
        else: 
            s.rho = rho_guess
            if s.rho.shape[0] != n: raise ValueError("Input HF guess has wrong shape")
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
            s.count += 1
            if s.count % 10 == 0: print "SCF Iteration:",s.count
            w,C = numpy.linalg.eigh(F)
            
            # Check for HOMO-LUMO degeneracy
            homo = w[s.m-1]
            lumo = w[s.m]
            if abs(lumo-homo) < 1e-4:
                kennyloggins.warning("DANGER ZONE: HOMO and LUMO are degenerate")

            Cocc = C[:,:s.m]
            rho_new = Cocc.dot(Cocc.T)

            ## Hubbard interaction
            if type(s.V) == float:
                F_new = s.h + numpy.diag(numpy.diag(rho_new)*s.V)
            elif len(s.V.shape)==4:
            ## General system; physicist's notation
                F_new = s.genFock(rho_new)
            else:
                print "Unknown V shape"
                import sys
                sys.exit(1)
            F_new_MO = C.T.dot(F_new).dot(C)
            #err = F_new_MO[:m,m:] #not sure what this would do instead?
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
#            emat = Cocc.T.dot(s.h+F).dot(Cocc) #for debugging at each iteration
#            kennyloggins.debug("Levels: %s" % str([emat[i,i] for i in range(len(emat))]))
        emat = Cocc.T.dot(s.h+F).dot(Cocc)
        E = 0
        for i in range(s.m):
            E += emat[i,i]
        s.F = F
        s.C = C
        s.Cocc = Cocc
        s.rho = Cocc.dot(Cocc.T)
        s.E = E

    def genFock(s,rho_new):
        F_new = s.h.copy()
        ri = range(len(s.h))
        for mu,nu,la,si in itertools.product(ri,ri,ri,ri):
            F_new[mu,nu] += (2*s.V[mu,si,nu,la] - s.V[mu,si,la,nu])*rho_new[la,si]
        return F_new

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

    def pdm2(s,P1):
        n = P1.shape[0]
        P2 = numpy.zeros((n,n,n,n))
        sz = range(n)
        for i,j,k,l in itertools.product(sz,sz,sz,sz):
            P2[i,j,k,l] = P1[k,i]*P1[j,l]
        return P2

    def mp2(s):
		E_mp2 = 0.
		for i in range(n_occ):
		 for j in range(n_occ):
		  for a in np.arange(n_occ,S_size,1):
		   for b in np.arange(n_occ,S_size,1):
			E_mp2 += MO_int[i,a,j,b]*(2*MO_int[i,a,j,b]-MO_int[i,b,j,a])/(eps[i]+eps[j]-eps[a]-eps[b])
		
		return E_mp2



if __name__ == "__main__":
    import est.lattice as lattice
    import est.estlib as estlib
    import sys
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', help='Number of electron pairs', type=int, required=True)
    parser.add_argument('-l', help='Lattice file', type=str, required=True)
#    parser.add_argument('-o', help='Output file', type=str, required=True)
#    parser.add_argument('-guess', help='HF Guess', type=str)
    args = parser.parse_args()

#    rho = numpy.genfromtxt(args.guess)

    lat = lattice.Lattice.fromfile(args.l,args.m)
    lat.V = estlib.EnforcePhysSym(lat.V)

#    if args.guess != None: fock = HF(lat.h,lat.V,lat.m,rho_guess = rho) 
#    else: fock = HF(lat.h,lat.V,lat.m)
    fock = HF(lat.h,lat.V,lat.m)
    fock._do_hf()
    print fock.E
#    numpy.savetxt(args.o,rho)


