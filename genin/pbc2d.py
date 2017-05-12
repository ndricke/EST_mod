import lattice
import numpy


import argparse

tol = 1e-8
#parser = argparse.ArgumentParser()
#parser.add_argument('-m', help="Number of electron pairs", type=int)
#parser.add_argument('-lat', help="Path to lattice file", type=str)
#args = parser.parse_args()
#m = args.m
#lat = args.lat

m_max = 648
lat_path = '2d/lats/'
#lat_list = ['x24y48u4p00.lat','x24y48u4p10.lat','x24y48u4p01.lat','x24y48u4p11.lat']
#lat_list = ['x36y36u4p11.lat']
lat_list = ['x36y36u4p00.lat','x36y36u4p10.lat','x36y36u4p01.lat','x36y36u4p11.lat']
print "m",' '.join(lat_list)

ldict = {} #dictionary mapping string to set of eigenvalues
for item in lat_list:
    l = lattice.Lattice.fromfile(lat_path+item,1)
    h = l.h
    e,v = numpy.linalg.eigh(h)
    ldict[item] = e


for m in range(1,m_max+1):
    hlgap_ttable = []
    for item in lat_list:
        e = ldict[item]
#        if m == 94: print e[m-3:m+3]

        homo_lumo_gap = e[m-1] - e[m]
        if abs(homo_lumo_gap) > tol:
            hlgap_ttable.append('1')
        else: hlgap_ttable.append('0')
    print m,' '.join(hlgap_ttable)



