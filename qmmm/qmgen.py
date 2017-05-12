import sys
import os
import gensphere
import inputlib
import numpy as np


def h2oTrim(pdb_in,outfile,cut_radius):
    output = []
    with open(pdb_in,'r') as f:
        for i,line in enumerate(f):
            output.append(line)
            if i == 4:
                c_coord = np.array([float(i) for i in line.split()[5:8]]) #determine central atom coord
                break
        for line in f:
            spline = line.split()
            if len(spline) > 4 and spline[2:4] == ["OW","SOL"]: #solvent water, measure from oxygen
                #r = float(spline[-1]) #supposedly traj order allows this, but Valerie mentioned it had errors
                coord = np.array([float(i) for i in spline[5:8]])
                r = np.linalg.norm(coord-c_coord)
                if r < cut_radius: #if the O is within the cutoff range
                    #print "Less than cutoff:",r
                    output.append(line)
                else: #need to delete the whole water molecule by skipping the following H lines
                    #print "Greater than cutoff:",r
                    f.next();f.next() #when I was previously reading the lines individually from file
            else:
                output.append(line)

    with open(outfile,'w') as f:
        f.writelines(output)

    return c_coord

def readxyz(filename):
    with open(filename,'r') as fin:
        data = fin.read().splitlines(True)[2:] #skip first 2 lines
    return data

def qchem2xyz(filename):
    f = open(filename)
    while True:
        line = f.readline()
        if not line:
            break
        if "Standard Nuclear Orientation" in line:
            f.readline()
            f.readline()
            lines = []
            while True:
                line = f.readline()
                if not line or "---------------------------------------" in line:
                    break
                lines.append(" ".join(line.split()[1:]) + '\n')
    return lines

    
class QmmmWrite():
    def __init__(self):
        self.conn = {}
        self.atom_dict = {}
        self.force_dict = {}
        self.mmIDdict = {'O':'186','H':'187','N':'23','C':'11','HN':'24','HC':'12'}

    def arShell(self,radius=16.0,N = 1000):
        self.Ar = N #need to keep track of this for figuring out how many atoms to freeze
        cen_atom = np.array([float(i) for i in self.qxyz[0].split()[1:4]])
        print cen_atom
        xc,yc,zc = gensphere.fiboSph(radius,cen_atom,N) #generate coordinates for Ar atoms around spherical shell
        self.ar_list = []
        for i in range(len(xc)):
            self.ar_list.append('Ar  %s  %s  %s  208   0 0 0 0' % (xc[i],yc[i],zc[i]))
        self.qxyz = self.qxyz + self.ar_list #need to do Ar second for how we index QM and waters

    #parse spline from pdb file. connect gives covalent connectivity preceded by force field ID
    def parseLinePdb(self,spline,mm_force,r_conn):#expects spline as list, connect as string
        atom = spline[2]
        if atom == "OW": atom = "O"
        if atom == "HW1" or atom == "HW2": atom = "H"
        coord = spline[5:8]
        app = atom+' '+' '.join(coord)+' '+mm_force+' '+r_conn 
        return app

    #read how atoms are connected from babel generated pdb file (before solvent added)
    #determines how atoms are connected(self.conn), and what MM force field for each atom (self.force_dict)
    def readConn(self,pdb_in):
        self.qm_count = 0 #counters number of QM atoms before we throw other junk in
        with open(pdb_in,'r') as f:
            for line in f:
                spline = line.split()
                if spline[0] == "HETATM" and spline[3] == "LIG": #check both, just to make sure we only take QM region
                    self.atom_dict[int(spline[1])] = spline[2]
                    self.qm_count += 1
                if spline[0] == "CONECT": #Yes, they spelled connect incorrectly
                    center = int(spline[1])
                    bonded = spline[2:] #all the atoms bonded to the central atom. We assume at least 1
                    if center not in self.conn: self.conn[center] = bonded
                    else: [self.conn[center].append(item) for item in bonded]
        #Now we determine what force field to use for each of the H's
        count = -1 #so that each user defined atom type is unique
        for key in self.atom_dict.keys():
            if self.atom_dict[key] == 'H':
                bond = self.atom_dict[int(self.conn[key][0])]  #conn[key][0] gives 1st entry of bonded atoms, a string
                if bond == 'N': self.force_dict[key] = '24'
                elif bond == 'C': self.force_dict[key] = '12'
            else: 
                try:
                    self.force_dict[key] = self.mmIDdict[self.atom_dict[key]] 
                except KeyError: #if the key isn't there, we make a new one
                    print "No entered force field for atom; resorting to user defined"
                    self.force_dict[key] = str(count)
                    count -= 1

#writes Q-chem xyz format from pdb data. As it reads the pdb, it uses the connectivity data stored earlier
    def pdb2qxyz(self,pdb_in):
        outlist = []
        count = 1 #count number of atoms so we don't need to check list repeatedly; starts at 1 in qchem
        with open(pdb_in,'r') as f:
            for line in f:
                spline = line.split()
                if len(spline) > 4:
                    if spline[3] == "LIG": #ligand molecule, treated quantum mechanically
                        try: 
                            conn_write = self.conn[count]
                            dl = len(conn_write) - 4 #Overlord Qchem demands connectivity to precisely 4 atoms
                            if dl > 0: conn_write = '0 0 0 0' #it's a metal, just leave it unbonded
                            elif dl < 0: conn_write.append('0 '*(-1*dl)) #add 0's to 4
                            else: pass #an atom with 4 bonds, which we need not change
                        except KeyError: #this probably means there is no connectivity, ie single atom
                            print "Connectivity Key Error; giving no connectivity"
                            conn_write = '0 0 0 0'
                        app = self.parseLinePdb(spline,self.force_dict[count],' '.join(conn_write))
                        outlist.append(app)
                        count += 1
                    elif spline[2:4] == ["OW","SOL"]: #solvent water, assumes order is O,H,H in file
                        con_O = ' '.join([str(count+1),str(count+2),'0 0'])
                        con_H = str(count)+' 0 0 0'
                        outlist.append(self.parseLinePdb(spline,'186',con_O))
                        outlist.append(self.parseLinePdb(f.next().split(),'187',con_H))
                        outlist.append(self.parseLinePdb(f.next().split(),'187',con_H))
                        count += 3
            assert count == len(outlist) + 1
        self.qxyz = outlist

    #Q-Chem QM/MM requires us to explicitly list every QM atom, which this function handles
    def qmWrite(self,f):
        f.write("\n$qm_atoms\n")
        f.write(' '.join([str(i) for i in range(1,self.qm_count+1,1)])+'\n') #write index for all QM atoms
        #write for self.qm_count+1 so that it doesn't cut off last atom
        f.write("$end\n\n")

#parse a filename for standardized charge and multiplicity information
#assumes last thing before .in is (a|c)(charge)(m)(multiplicity)
    def chmult(self,fname):
        spl = fname.split('_')[-1]
        spl = spl.split('.')[0]
        assert len(spl) == 4 #this should be true unless charge or mult is over 10, or got wrong thing
        charge = int(spl[1])
        if spl[0] == 'a': charge *= -1
        mult = int(spl[-1])
        return (charge,mult)

    def wInpchmult(self,f,charge,mult):
        f.write("model_system_mult   " + str(mult)+'\n')
        f.write("model_system_charge " + str(charge)+'\n')

#Converts Q-chem QM/MM velocity file (scratch:AIMD/NucVeloc) for input; for use in wrapper script
    def qmVel(self,f):
        vel = np.genfromtxt(f,skip_header=1) #first line just says what the data is
        vel = vel[0,:][1:] #I think the first entry is a time step, followed by actual velocity data
        vel = vel.reshape((vel.shape[0]/3,3))
        for i in self.fix_list: #freeze QM atoms
            vel[i-1,:] = 0.0 #fix_list starts from 1
        self.vel = vel

    def writeVel(self,f):
        f.write("\n$velocity\n")
        for i in range(self.vel.shape[0]):
            f.write("%s %s %s\n" % (self.vel[i,0],self.vel[i,1],self.vel[i,2]))
        f.write("$end\n")

#Write out the list of fixed atoms. Should be QM atoms and Ar shell
    def qmFix(self,f):
        f.write("\n$opt\nFIXED\n")
        for i in self.fix_list:
            f.write(str(i) + " XYZ\n")
        f.write("ENDFIXED\n$end\n")

    #makes the qmmm input file
    def QmmmQin(self,outfile,vel_file):
        charge,mult = self.chmult(outfile)
        self.xyz_count = len(self.qxyz)
        self.qm_h2o_count = self.xyz_count - self.Ar #difference of total and Ar atoms is the qm+h2o atoms
        self.fix_list = [i for i in range(1,self.qm_count+1,1)]
        [self.fix_list.append(i) for i in range(self.qm_h2o_count+1,self.xyz_count+1,1)]
        self.qmVel(vel_file)

        with open(outfile,'w') as f:
            f.write(inputlib.qmmm_rem1)
            self.wInpchmult(f,charge,mult)
            f.write("aimd_fixed_atoms  " + str(self.qm_count+self.Ar)) #do we actually want to fix qm atoms?
            f.write("\n$end\n")
#            f.write(inputlib.forceman)
            f.write(inputlib.Fe_ff)
            self.qmWrite(f)
            f.write("$molecule\n")
            f.write(' '.join([str(charge),str(mult)])+"\n")
            f.write('\n'.join(self.qxyz)+'\n')
            f.write("$end\n\n")
            self.writeVel(f)
            self.qmFix(f)

            f.write(inputlib.read) #this includes the @@@ for between jobs
            f.write(inputlib.qmmm_rem2)
            self.wInpchmult(f,charge,mult)
            f.write("aimd_fixed_atoms  " + str(self.qm_count+self.Ar))
            f.write("\n$end\n")
#            f.write(inputlib.forceman)
            f.write(inputlib.Fe_ff)
            self.qmWrite(f)
            self.qmFix(f)
            self.writeVel(f)
        
    def qmVelGen(self, outfile):
        charge,mult = self.chmult(outfile)
        self.xyz_count = len(self.qxyz)
        self.qm_h2o_count = self.xyz_count - self.Ar #difference of total and Ar atoms is the qm+h2o atoms
        self.fix_list = [i for i in range(1,self.qm_count+1,1)]
        [self.fix_list.append(i) for i in range(self.qm_h2o_count+1,self.xyz_count+1,1)]
        with open(outfile,'w') as f:
            f.write(inputlib.qmmm_remvel)
            self.wInpchmult(f,charge,mult)
            f.write("aimd_fixed_atoms  " + str(self.qm_count+self.Ar)) #do we want to fix qm atoms?
            f.write("\n$end\n")
#            f.write(inputlib.forceman) #I think this is for pbc
            f.write(inputlib.Fe_ff)
            self.qmWrite(f)
            f.write("$molecule\n")
            f.write(' '.join([str(charge),str(mult)])+"\n")
            f.write('\n'.join(self.qxyz)+'\n')
            f.write("$end\n\n")
            self.qmFix(f)

    def frzIn(self,outfile):
        self.charge,self.mult = self.chmult(outfile)
        self.xyz_count = len(self.qxyz)
        self.qm_h2o_count = self.xyz_count - self.Ar #difference of total and Ar atoms is the qm+h2o atoms
        self.fix_list = [i for i in range(1,self.qm_count+1,1)]
        [self.fix_list.append(i) for i in range(self.qm_h2o_count+1,self.xyz_count+1,1)]

        with open(outfile,'w') as f:
            self.remWrite(f,inputlib.qmmm_remvel)
            self.randWrite(f)
            f.write("$molecule\n")
            f.write(' '.join([str(self.charge),str(self.mult)])+"\n")
            f.write('\n'.join(self.qxyz)+'\n')
            f.write("$end\n\n")

            f.write(inputlib.read) #this includes the @@@ for between jobs
            self.randWrite(f)
            self.remWrite(f,inputlib.qmmm_rem1)

            f.write(inputlib.read) #this includes the @@@ for between jobs
            self.randWrite(f)
            self.remWrite(f,inputlib.qmmm_rem2)
   
    def remWrite(self,f,rem):
        f.write(rem)
        self.wInpchmult(f,self.charge,self.mult)
        f.write("aimd_fixed_atoms  " + str(self.qm_count+self.Ar))
        f.write("\n$end\n")
    
    def randWrite(self,f):
        f.write(inputlib.Fe_ff)
        self.qmWrite(f)
        self.qmFix(f)

    def xyzWrite(self):
        with open('q.xyz','w') as f:
            f.write(str(len(qmmm.qxyz))+'\n')
            f.write("Energy:    100"+'\n')
            for line in qmmm.qxyz:
                spline = line.split()[:4]
                f.write(' '.join(spline)+'\n')







