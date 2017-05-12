import numpy
import itertools
import xform
import sys

class DMETEnergy:
    def _calc_energy(s,frag,fock):
        P1f = frag.solver.get_P() #fragment 1pdm
        center = [int(i) for i in frag.center] #convert center names to fragment index, rather than lattice index
        frag.e1 = 0
        frag.e1c = 0

        hf = frag.h
        nz = range(2*frag.ni)
        niz = range(frag.ni) #we always organize fragment, then bath

        for f,j in itertools.product(niz,nz):
            frag.e1 += hf[f,j] * P1f[j,f]
        for f,j in itertools.product(center,nz):
            frag.e1c += hf[f,j] * P1f[j,f]

#        #Old way of calculating core energy contribution, need to also modify bootstrap part
#        frag.core = 0.0; frag.core_cen = 0.0
#        for f in niz:
#            frag.core += frag.h[f,f]*P1f[f,f]
#        for c in center:
#            frag.core_cen += frag.h[c,c]*P1f[c,c]
#        for j in nz:
#            for f in niz:
#                if f != j:
#                    frag.e1 += hf[f,j] * P1f[j,f]
#            for f in center:
#                if f != j:
#                    frag.e1c += hf[f,j] * P1f[j,f]
#        frag.e1 += 0.5*frag.core
#        frag.e1c += 0.5*frag.core_cen

        frag.e1 /= frag.ni
        frag.e1c /= len(center)

        #2 electron component
        if type(s.lat.V) == float: frag.e2,frag.e2c = s.e2_Ucalc(frag)
        elif len(s.lat.V.shape) == 4: frag.e2,frag.e2c = s.e2_V4calc(frag,center,fock)
        frag.e = frag.e1 + frag.e2 
        frag.ec = frag.e1c + frag.e2c 

#Calculates 2e- energy when V is a 4-tensor
    def e2_V4calc(s,frag, center,fock):
        e2 = 0.0
        e2c = 0.0
        f_index = frag.lattice_indices() 
        niz = range(2*frag.ni) 
        P2f = frag.solver.get_P2()
        
        #Fragment 2e Energy
        for f,j,k,l in itertools.product(f_index,niz,niz,niz):
            e2 += frag.V[f,j,k,l]*P2f[f,j,k,l]
        e2 /= frag.ni*4.
        
        #Center 2e Energy
        for c,j,k,l in itertools.product(center,niz,niz,niz):
            e2c += frag.V[c,j,k,l]*P2f[c,j,k,l]
        e2c /= len(center)*4.

        return e2,e2c

#Calculates 2e- energy when V is a hubbard U, so we assume V is diagonal
    def e2_Ucalc(s,frag):
        center_fragdex = [int(i) for i in frag.center] #we're working within the fragment numbering
        P2_s = frag.solver.get_P2() #as this 2pdm is only in the Schmidt space
        e2 = 0.0; e2c = 0.0
        for i in range(frag.ni):
            e2 += P2_s[i,i,i,i]/4.*frag.V[i,i,i,i]
        for i in center_fragdex:
            e2c += P2_s[i,i,i,i]/4.*frag.V[i,i,i,i]
        e2 /= frag.ni
        e2c /= len(center_fragdex)
        return e2,e2c

    def energy(s):
        return s.e

    def center_energy(s):
        return s.ec

    def pr_P2(s,diag=True,off_diag=False):
        if diag == True:
            print "P2 Diagonal Entries:"
            for i in range(s.ni):
                print s.P2_r[i,i,i,i]
        if off_diag == True:
            print "P2 Off-Diagonal Entries:"
            for i in range(s.ni-1):
                print str(s.P2_r[i,i+1,i,i+1]) + " " + str(s.P2_r[i+1,i,i+1,i]) 



