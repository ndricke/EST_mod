import ChemData as CD
import sys
import numpy as np

#Single point/general job parse keys
finish_key ='Thank you very much for using Q-Chem.  Have a nice day.'
sp_key ='Total energy in the final basis set'
coord_key ='Standard Nuclear Orientation (Angstroms)'

#Freq job keys
freq_key ='Frequency:   '
entropy_key ='Total Entropy:'
enthalpy_key ='Total Enthalpy:'
zpe_key = 'Zero point vibrational energy:'
vibH_key ='Vibrational Enthalpy:'

#PCM job keys
gelec_key ='G_electrostatic'
solG_key ='Total Free Energy (H0 + V/2 + non-elec)'
solEtot_key ='Total energy in the final basis set'
job_key = 'jobtype'
chmul_key = '$molecule'

#Parsing class for Q-Chem output files
class Qdata(object):
  def __init__(self):
    #There was a point where I didn't include all of these, but hey, why not? It runs fast enough
    self.parse_base = {finish_key:self.isComplete, \
                       chmul_key:self.chargeMult, \
                       sp_key:self.spEnergy, \
                       coord_key:self.coordGrab, \
                       entropy_key:self.entropy,\
                       enthalpy_key:self.enthalpy,\
                       zpe_key:self.zeroPointEnergy,\
                       vibH_key:self.enthalpyVib, \
                       gelec_key:self.solvationG,
                       solG_key:self.solvEnergy, \
                       job_key:self.jobtype, \
                       }

#main workhorse; iterate through file and grab stuff when it hits a key phrase
  def qParse(self, qchem_outfile):
    with open(qchem_outfile, 'r') as f:
      for line in f:
        #if match with parse list, call function
        for item in self.parse_base.keys():
          if item in line:
            func_point = self.parse_base[item]
            func_point(f,line)

  def readFile(self, filename):
    spl_filename = filename.split('.')
    if spl_filename[-1] == 'xyz': #Just get the coordinates from the xyz file
      with open(filename,'r') as f:
        coord_list = f.read().splitlines(True)[2:]
        self.coordArr(coord_list)
    elif spl_filename[-1] == 'out': #Q-Chem output file
      self.qParse(filename)

  def default(self):
    self.basis = '6-31+g*'
    self.method = 'b3lyp'
    self.mult = 1
    self.charge = 0
    self.job = 'sp'

  def listCoord(self):
    coord_list = []
    for i,atom in enumerate(self.atoms):
      coord_list.append(' '.join([atom] + [str(co) for co in list(self.coord[i,:])]))
    return coord_list

  def isComplete(self, infile = None,line=None):
    self.finish = True

  def chargeMult(self,infile,line):
    chmul = infile.next().strip('\n').split()
    self.charge = int(chmul[0])
    self.mult = int(chmul[1])
    del self.parse_base[chmul_key]

  def spEnergy(self, infile, line):
    spline = line.split()
    self.E = float(spline[-1])

  def solvEnergy(self,infile,line):
    spline = line.split()
    self.E = float(spline[-2])

  def entropy(self, infile, line):
    spline = line.split()
    self.S = float(spline[2])/1000.0
    del self.parse_base[entropy_key]

  def enthalpy(self, infile, line):
    spline = line.split()
    self.H = float(spline[2])
    del self.parse_base[enthalpy_key]

  def enthalpyVib(self, infile, line):
    spline = line.split()
    self.Hvib = float(spline[2])
    del self.parse_base[vibH_key]

  def frequencies(self, infile, line):
    spline = line.split()
    [self.freq.append(float(frq)) for frq in spline[1:]]

  def zeroPointEnergy(self, infile, line):
    spline = line.split()
    self.Hzpe = float(spline[-2])
    del self.parse_base[zpe_key]

  def solvationG(self, infile, line):
    spline = line.split()
    self.Gpcm = float(spline[-2])
    self.solvation = 'pcm'

  def coordGrab(self, infile, line):
    key_end = "--------------"
    infile.next() ; infile.next() #Skipping junk. 2nd line contains ------
    grab_list = []
    for line in infile:
      if key_end in line:
        break
      grab_list.append(line.strip())
    self.coordArr(grab_list)
  
  def jobtype(self, infile, line):
    spline = line.split()
    self.job = spline[-1]

  def coordArr(self,coord_list):
    self.coord = np.zeros((len(coord_list),3))
    self.atoms = []
    for i, line in enumerate(coord_list):
      spline = line.split()
      self.atoms.append(spline[1])
      self.coord[i,:] = [float(j) for j in spline[2:]]

    

#short testing module for reading a frequency calculation
if __name__ == "__main__":
  qdata = Qdata(sys.argv[1])
  print qdata.freq


















