import numpy as np
import itertools
import copy

import schmidt
import xform
import hf
import dmetenergy

np.set_printoptions(precision=4,suppress=True)

#Numerical recipes implementation of bracket method to find sufficient range for finding a root
def bracket(f,xbrack,maxiter=50,xconst=1.6):
    if xbrack[0] == xbrack[1]: raise ValueError("guess should be 2 different values")
    f0 = f(xbrack[0]); f1 = f(xbrack[1]);
    for i in range(maxiter):
        if f0*f1 < 0: return xbrack,True
        elif abs(f0) < abs(f1):
            xbrack[0] += xconst*(xbrack[0]-xbrack[1])
        else:
            xbrack[1] += xconst*(xbrack[1]-xbrack[0])
    return xbrack,False

def antisymMat(n):
    Z = (np.random.rand(n,n)-0.5)*0.5
    Z = Z+Z.T
    for i in range(n):
        Z[i,i] = 0.
        for j in range(i):
            Z[i,j] *= -1.
    return Z

def unitMatGen(n):
    Z = antisymMat(n)
    U = scl.expm(Z)
    B = 2.*(np.random.rand(n)-0.5)
    B = np.diag(B)
    A = U.dot(B.dot(U.T))
    return Z,U,B,A

def redFock(n,rho,h,V):
    F_red = h.copy()
    ri = range(n)
    for mu,nu,la,si in itertools.product(ri,ri,ri,ri):
        F_red[mu,nu] += (V[mu,si,nu,la] - 0.5*V[mu,si,la,nu])*rho[la,si]
#    F_red = c1c2.hRedHalf
    return F_red

def embed_energy(nf,h,V,rho): 
    n = len(h) 
    ri = range(n) 
    rf = range(nf) 
     
    e1 = 0 
    for f,j in itertools.product(rf,ri): 
        e1 += 2.*h[f,j]*rho[j,f] 
 
    e2 = 0 
    e2 = engV(n,nf,V,rho)
    print "e1: ", e1
    print "e2: ", e2

    return e1 + e2

def engV(n,nf,V,rho):
    ri = range(n) 
    rf = range(nf) 
    e2 = 0.
    for f,j,k,l in itertools.product(rf,ri,ri,ri): 
        e2 += (2*V[f,j,k,l]-V[f,j,l,k])*rho[f,k]*rho[j,l] 
    return e2

def hf_energy(h,V,rho):
    e1 = np.trace(h.dot(rho)) * 2
    e2 = 0

    n = len(h)
    ri = range(n)
    gamma = np.zeros((n,n,n,n))
    for i,j,k,l in itertools.product(ri,ri,ri,ri):
        gamma[i,j,k,l] = rho[k,i]*rho[l,j]
    for i,j,k,l in itertools.product(ri,ri,ri,ri):
        e2 += (2*V[i,j,k,l]-V[i,j,l,k])*gamma[k,l,i,j]

    print "e1: ", e1
    print "e2: ", e2
    return e1 + e2

def Vhubbard(n, U):
    V = np.zeros((n,n,n,n))
    for i in range(n):
        V[i,i,i,i] = U
    return V

def hRedHalf(nf,h,T):
    h_r = copy.copy(h)
    h_r[nf:,:] = 0.
    return xform.one(T,h_r)

def VRedHalf(nf,V,T):
    V_r = copy.copy(V)
    V_r[nf:,:,:,:] = 0.
    return xform.two(T,V_r)

def hRedSym(nf,h,T):
    h_r = copy.copy(h)
    h_r[nf:,:nf] *= 0.5
    h_r[:nf,nf:] *= 0.5
    h_r[nf:,nf:] *= 0.0
    return xform.one(T,h_r)

def hubb(n,U):
    V = Vhubbard(n,U)
    h = np.diag(np.ones(n-1),1)
    h[0,-1] = -1
    h += h.T
    h *= -1
    return h,V

def randHV(n,U=1.,scale=0.1):
    h,V = hubb(n,U)
    h_rand = (np.random.rand(n,n)-0.5)*scale
    h_rand += h_rand.T
    h += h_rand
    V_rand = np.zeros((n,n,n,n))
    nz = range(n)
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        V_rand[i,j,k,l] = np.random.rand()*scale*U
        V_rand[j,i,l,k] = V_rand[i,j,k,l]
        V_rand[k,l,i,j] = V_rand[i,j,k,l]
        V_rand[l,k,j,i] = V_rand[i,j,k,l]
        V_rand[k,j,i,l] = V_rand[i,j,k,l]
        V_rand[l,i,j,k] = V_rand[i,j,k,l]
        V_rand[i,l,k,j] = V_rand[i,j,k,l]
        V_rand[j,k,l,i] = V_rand[i,j,k,l]
    V += V_rand
    return h,V

def C1():

    return C1

def C2():
    
    return C2

if __name__ == "__main__":
    n = 12
    m = 4
    U = 2.
    nf = 2

    h,V = randHV(n,U=U,scale=0.)

    print "Beginning HF"
    hart = hf.HF(h,V,m)
    psi = hart.get_Cocc()
    rho = hart.get_rho()
    print hart.get_energy()/n

    #Method 1, producing C1
    M = psi[:nf,:] #select the first nf orbitals out of the HF wavefunction
    S_M = np.dot(M.T,M) #calculate the overlap matrix of the fragment orbitals
    u,s,v = np.linalg.svd(M) #perform an SVD of the fragment overlap matrix
    C = psi.dot(v.T) #wavefunction dot SVD
    T = np.zeros((n,2*nf)) #pull components out of C; split entangled bath and fragment into block diagonal matrix
    T[:nf,:nf] = C[:nf,:nf]
    T[nf:,nf:] = C[nf:,:nf]
    Tc = C[:,nf:]

    C1 = copy.copy(T)
    #C1[:nf,:nf] = np.eye(nf)

    for i in range(nf):
        C1[:,nf+i] /= np.sqrt(C1[:,nf+i].dot(C1[:,nf+i]))
        C1[:,i] /= np.sqrt(C1[:,i].dot(C1[:,i]))

    h_f1 = xform.one(C1,h)
    V_f1 = xform.two(C1,V)
    core = xform.core(Tc,V)
    rho_f1 = xform.one(C1,rho+core)

    print embed_energy(nf,h_f1,V_f1,rho_f1)/nf

    #Method 2 for producing C2, which is a unitary rotation of C1
    rho_f = rho[:nf,:nf]
    Gam,U = np.linalg.eigh(rho_f)
    A = U.dot(np.diag(Gam**0.5))
    ATinv = np.linalg.inv(A.T)
    B = rho[nf:,:nf].dot(ATinv)

    rhoh = np.eye(n) - rho
    rhoh_f = rhoh[:nf,:nf]
    Gamh,Uh = np.linalg.eigh(rhoh_f)
    Ah = Uh.dot(np.diag(Gamh**0.5))
    AhTinv = np.linalg.inv(Ah.T)
    Bh = rhoh[nf:,:nf].dot(AhTinv)

    C2 = np.zeros((n,2*nf))
    C2[:nf,:nf] = A
    C2[:nf,nf:] = Ah
    C2[nf:,:nf] = B
    C2[nf:,nf:] = Bh

    h_f2 = xform.one(C2,h)
    V_f2 = xform.two(C2,V)
    rho_f2 = xform.one(C2,rho+core)

    h_r0 = hRedHalf(nf,h,C2)
    V_f2_0 = VRedHalf(nf,V,C2)

    h_r05 = hRedSym(nf,h,C2)

    print (2.*np.trace(h_r0.dot(rho_f2)) + engV(2*nf,2*nf,V_f2_0,rho_f2))/nf
    print (2.*np.trace(h_r05.dot(rho_f2)) + engV(2*nf,2*nf,V_f2_0,rho_f2))/nf
















