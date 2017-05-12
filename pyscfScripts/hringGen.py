import numpy as np
import sys

#This class generates the vertices of an n-sided polygon inscribed in a circle.
#The radius is implicitly defined by side_length, which is the distance between vertices
class InscPoly(object):
    def __init__(self,polygon_n,side_length):
        self.n = polygon_n
        self.a = side_length
        self.r = self.radius()
        self.geom = self.coords()

    def radius(self):
        return self.a/(2.0*np.sin(np.pi/self.n))

#generate points for the (presently regular) polygon, center at origin
    def coords(self):
        rad_list = np.arange(self.n)*2.0*np.pi/self.n
        coord_list = []
        for rad in rad_list:
            x = self.r*np.cos(rad)
            y = self.r*np.sin(rad)
            coord_list.append((str(x),str(y),'0.0'))
        return coord_list

#This inherits from the InscPoly class because an inscribed polygon is the location of the H's in an H-ring
class Hring(InscPoly):
    def __init__(self,polygon_n,side_length):
        super(Hring, self).__init__(polygon_n,side_length) #generate coordinates for H-ring
        #Set two variables for generating Q-Chem input files
        self.molecule = \ 
"""$molecule
0 1
"""
        self.rem = \
"""$rem
method hf
basis sto-3g
correlation 105
scf_guess core
symmetry false
sym_ignore true
purecart 111
iprint 200
print_orbitals=99999
use_new_path2 false
save_ao_integrals 2
$end
"""
    
    #Writes Q-Chem input file for H-ring
    def qchemHring(self, qinfile):
        with open(qinfile,'w') as f:
            f.write(self.molecule)
            for item in self.geom:
                f.write("H "+" ".join(item)+"\n")
            f.write("$end\n\n")
            f.write(self.rem)


if __name__ == "__main__":

    #Hbond = 0.74 #distance for H2 bond
    Hbond = 0.951 #distance for 10H ring at equilibrium

    n = int(sys.argv[1]) #number of H's in H-ring
    infilename = sys.argv[2] #name of Q-Chem input file that is generated

    poly = Hring(n,Hbond) #construct H-ring with n H's at distance Hbond from each other
    poly.qchemHring(infilename) #write a Q-Chem input file for solving HF and extracting integrals

