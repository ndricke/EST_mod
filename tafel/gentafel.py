import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


#prev = partially reversible cycle, where every step is reversible except for one
#this function assumes that step 4 is irreversible
def prevCycle4(k1,k2,k3,k4,kn1,kn2,kn3):
    num = k1*k2*k3
    denom = kn1*kn2*kn3+ kn1*kn2*k4+ kn1*k3*k4+ k2*k3*k4+ kn2*kn3*k1+ kn2*k1*k4+ k1*k3*k4+ \
            kn3*k1*k2+ k1*k2*k4+ k1*k2*k3
    c4 = num/denom
    r4 = k4*c4
    c3 = c4*(kn3+k4)/k3
    c2 = (kn2*c3+k4*c4)/k2
    c1 = (kn1*c2+k4*c4)/k1
    return r4, [c1,c2,c3,c4]

def prevCycle2(k1,k2,kn1):
    c2 = k1/(k1+kn1+k2)
    r2 = k2*c2
    c1 = c2*(kn1+k2)/k1
    return r2, [c1,c2]




#should be given 2n-1 rate constants, where n is the number of steps in the cycle
def genPrevCycle(k):
    pass

def prevOH(OH,ks,OHdep,prevFunc):
    for i in range(len(ks)):
        ks[i] *= OH**OHdep[i]
    return prevFunc(*ks)


font = {'size':18}
mpl.rc('font',**font)

#Want to scan over pOH dependent steps for log(j) dependence on current
#Let's assume that all of the reversible steps are pOH dependent in the reverse direction

dpoints = 500
OH = np.linspace(0.01,5,dpoints)
OHdep = np.array([0,0,0,0,1,1,1])
ks = np.array([1.,1.,1.,1.,1.,1.,1.])

prev = prevCycle4
r_arr = np.zeros(dpoints)
c = np.zeros((dpoints,4))
for i in range(dpoints):
    r, c_list = prevOH(OH[i],ks,OHdep,prev)
    r_arr[i] = r
    c[i,:] = c_list

#OHdep = [0,0,1]
#ks = [1.,1.,1.]
#prev = prevCycle2
#r, c_list = prevOH(OH,ks,OHdep,prev)

#print c_list
for i in range(4):
    plt.plot(OH,c[:,i],label=r'$\theta_%s$ ' % str(i))
plt.xlabel('[OH]')
plt.ylabel('Surface Coverage')
plt.legend(loc=2)
plt.title("All k's = 1")
plt.show()

#print r
#print c_list


