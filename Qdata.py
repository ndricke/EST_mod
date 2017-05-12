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

#Parsing class for Q-Chem output files
class Qdata(object):
  def __init__(self, outfile):
    self.coord = None
    self.atoms = []
    self.E = None
    self.freq = []
    self.Gpcm = None
    self.finish = False
    self.H = None
    self.Hzpe = None
    self.parse_base = {finish_key:self.isComplete, \
                       sp_key:self.spEnergy, \
                       coord_key:self.coordGrab
                       }
    self.qParse(outfile)

#main workhorse; iterate through file and grab stuff when it hits a key phrase
  def qParse(self, qchem_outfile):
    with open(qchem_outfile, 'r') as f:
      for line in f:
        #if match with parse list, call function
        for item in self.parse_base.keys():
          if item in line:
            func_point = self.parse_base[item]
            func_point(f,line)

##I had thought I might use this code to automatically decide what to parse, but we do that from filename
##through gform instead
#  def inpParam(self, infile):
#    job = None
#    jobtype = infile.split('_')[-2]
#    inp_list = CD.sectionGrab('$molecule','-----------',infile, 0,1)

  def loadFq(self,fq_file):
    self.parse_base = {
                #freq_key:self.frequencies,\
                entropy_key:self.entropy,\
                enthalpy_key:self.enthalpy,\
                zpe_key:self.zeroPointEnergy,\
                vibH_key:self.enthalpyVib
               }
    self.qParse(fq_file)

  def loadPCM(self, pcm_file):
    self.parse_base = {gelec_key:self.solvationG, solG_key:self.solvEnergy, solEtot_key:self.spEnergy}
    self.qParse(pcm_file)

#Load information from another Qdata instance
  def readQdat(self,r_qdat,dtype):
    if dtype == 'fq':
      self.H = r_qdat.H
      self.S = r_qdat.S
      self.Hvib = r_qdat.Hvib
      self.Hzpe = r_qdat.Hzpe
    if dtype == 'sp':
      self.E = r_qdat.E

  def isComplete(self, infile = None,line=None):
    self.finish = True

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

  def coordGrab(self, infile, line):
    key_end = "--------------"
    infile.next() ; infile.next() #Skipping junk. 2nd line contains ------
    grab_list = []
    for line in infile:
      if key_end in line:
        break
      grab_list.append(line.strip())
    if self.coord == None: self.coord = np.zeros((len(grab_list),3))
    self.atoms = []
    for i, line in enumerate(grab_list):
      spline = line.split()
      self.atoms.append(spline[1])
      self.coord[i,:] = [float(j) for j in spline[2:]]

#short testing module for reading a frequency calculation
if __name__ == "__main__":
  qdata = Qdata(sys.argv[1])
  print qdata.freq


















