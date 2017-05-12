import numpy as np
import est.hf as hf
import est.xform as xform
import est.estlib as estlib


def qChemInt(indir):
    ao_ind = indir+"/801.0.txt"
    ao_2ei = indir+"/21.0.txt"
    mo_file = indir+"/53.0.txt"

    outfile = indir.split('/')[-1]
    qout = indir+'/'+outfile+'.out'

    #read indices and integrals into numpy arrays w/ sparse format
    inds = np.genfromtxt(ao_ind,skip_header=1)
    inds = inds.astype(int)
    ei2 = np.genfromtxt(ao_2ei)
    C = np.genfromtxt(mo_file)

    #turn to dense 4-tensor with physicist's notation (Qchem starts with chemist's notation)
    V = estlib.toDense(inds,ei2,"chem")
    V = estlib.chem2phys(V)
    n = V.shape[0]

    #Parse the qchem output file for S and h
    S_key = "Overlap Matrix"
    h_key = "Core Hamiltonian Matrix"

    with open(qout,'r') as f:
        for line in f:
            if S_key in line:
                S = estlib.matrixGrab(f,n)
            if h_key in line:
                h = estlib.matrixGrab(f,n)

    #np.savetxt('S.txt',S) #just in case we want it later

    #Symmetric orthogonalization of h and V for use in DMET code
    eS, vS = np.linalg.eigh(S)
    eS_nroot = np.diag(eS**-0.5)
    S_nroot = vS.dot(np.dot(eS_nroot,vS.T))

    eS_root = np.diag(eS**0.5)
    S_root = vS.dot(np.dot(eS_root,vS.T))

    h_o = xform.one(S_nroot,h)
    V_o = xform.two(S_nroot,V)

    ##Enforce symmetry of 1 and 2 electron integrals
    V_os = estlib.EnforcePhysSym(V_o)
    h_os = (h_o+h_o.T)/2.0 

    n = int(h_o.shape[0])
    m = int(n//2)

    ##Reshape the MO coefficients, and form the density matrix
    C = C[:(n*n)]
    C = C.reshape(n,n).T
    Cocc = C[:,:m]
    C_o = np.dot(S_root,Cocc)
    rho_o = np.dot(C_o,C_o.T)

    ##Code for testing whether integrals have been read properly
    fock = hf.HF(h_os,V_os,m,rho_guess=rho_o)
    F = fock.genFock(rho_o)
    print "1e-: ",np.trace(np.dot(h_os,rho_o))
    print "Total: ",np.trace(np.dot(h_os+F,rho_o))

    return n,h_os,V_os,m,rho_o








