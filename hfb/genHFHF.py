import numpy as np

import hf
import xform
import schmidt
import c1c2





n = 12
m = 4
U = 2.
nf = 3

#make a list of range(nf) that tile n
fraglist = []
for i in range(n/nf): #assumes n divisible by nf
    fraglist.append(range(nf*i,nf*(1+i)))
print fraglist

h,V = c1c2.randHV(n,U=U,scale=0.5)

print "Beginning HF"
hart = hf.HF(h,V,m)
psi = hart.get_Cocc()
rho = hart.get_rho()
print hart.get_energy()/n

e_list = []
for frag_sites in fraglist:
    T = schmidt.schmidt(psi,frag_sites,allv=True)
    C1 = T[:,:2*nf]
    Tc = T[:,2*nf:]
    h_f1 = xform.one(C1,h)
    V_f1 = xform.two(C1,V)
    core = xform.core(Tc,V)
    rho_f1 = xform.one(C1,rho+core)

    efrag = c1c2.embed_energy(nf,h_f1,V_f1,rho_f1)/nf
    e_list.append(efrag)
    print "Efrag %s: %s" % (" ".join([str(item) for item in frag_sites]),str(efrag))

e_arr = np.array(e_list)
print np.sum(e_arr)/(n/nf)



