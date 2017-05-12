import numpy as np
import itertools

def traceEnergy(h,V,P1,P2):
    e1 = 0.0; e2 = 0.0
    n = h.shape[0]
    nz = range(n)
    for i,j in itertools.product(nz,nz):
        e1 += h[i,j] * P1[j,i]
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        e2 += V[i,j,k,l] * P2[i,j,k,l]
    return e1+e2/4, e1, e2/4

def traceRG(h,V,P1,P2):
    e1 = 0.0; e2 = 0.0
    n = h.shape[0]
    nz = range(n)
    for i,j in itertools.product(nz,nz):
        e1 += h[i,j] * P1[j,i]
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        e2 += V[i,j,k,l] * P2[i,j,l,k]
    return e1+e2, e1, e2

def P2mfTrace(P2,V,ni):
    e2 = 0
    n = V.shape[0]
    nz = range(n); niz = range(ni)
    for i,j,k,l in itertools.product(niz,nz,nz,nz):
        e2 += (2.0*V[i,j,k,l]-V[i,j,l,k])*P2[k,l,i,j]
    return e2

def chem2phys(P2c):
    n = P2.shape[0]
    P2p = np.zeros((n,n,n,n))
    nz = range(n)
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        P2p[i,j,k,l] = P2c[i,k,j,l]
    return P2p

def P2troy2gen(P2troy,V):
    thresh = 10**-8
    n = P2troy.shape[0]
    P2gen = np.zeros((n,n,n,n))
    nz = range(n)
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        denom = (2*V[i,j,k,l] - V[i,j,l,k])
        if denom <= thresh: P2gen[k,l,i,j] = 0.
        else:
            P2gen[k,l,i,j] = V[i,j,k,l]/denom*P2troy[i,j,k,l]
    return P2gen

def P2gen2troy(P2gen,V):
    thresh = 10**-8
    n = P2gen.shape[0]
    P2troy = np.zeros((n,n,n,n))
    nz = range(n)
    for i,j,k,l in itertools.product(nz,nz,nz,nz):
        denom = V[i,j,k,l]
        if denom <= thresh: P2troy[k,l,i,j] = 0.
        else:
            P2troy[i,j,k,l] = (2*V[i,j,k,l]-V[i,j,l,k])/denom*P2gen[i,j,k,l]
    return P2troy

def read2Tensor(tensor_file):
    with open(tensor_file,'r') as f:
        n = int(f.readline())
        tensor = np.zeros((n,n))
        for line in f:
            s = line.split()
            i = int(s[0]); j = int(s[1]);
            v = float(s[2])
            tensor[i,j] = v
    return tensor

def read4Tensor(tensor_file):
    with open(tensor_file,'r') as f:
        n = int(f.readline())
        tensor = np.zeros((n,n,n,n))
        for line in f:
            s = line.split()
            i = int(s[0]); j = int(s[1]); k = int(s[2]); l = int(s[3]);
            v = float(s[4])
            tensor[i,j,k,l] = v
    return tensor

def calcEdmet(h,V,P1,P2,P2_mf,ni):
    e1 = 0; e2 = 0; e2_fb = 0; e2_mf = 0
    niz = range(ni)
    nz = range(h.shape[0])
    for f,j in itertools.product(niz,nz):
        e1 += h[f,j] * P1[j,f]
    e1 /= ni

    for i,j,k,l in itertools.product(niz,nz,nz,nz):
        e2_fb += V[i,j,k,l]*P2[i,j,k,l]
        e2_mf += (2.0*V[i,j,k,l]-V[i,j,l,k])*P2_mf[k,l,i,j]
    e2_fb /= ni*4.
    e2_mf /= ni*2.
    e2 = e2_fb + e2_mf
    return e1+e2,e1,e2,e2_mf

def dmetEfix(h,V,P1,P2,ni):
    e1 = 0; e2 = 0;
    niz = range(ni)
    nz = range(h.shape[0])
    for f,j in itertools.product(niz,nz):
        e1 += h[f,j] * P1[j,f]
    for i,j,k,l in itertools.product(niz,nz,nz,nz):
        e2 += V[i,j,k,l]*P2[i,j,k,l]
    e1 /= ni
    e2 /= ni
    return e1+e2,e1,e2





