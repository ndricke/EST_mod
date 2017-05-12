import Qdata
import ChemData
import pandas as pd
import sys
import os
from collections import Counter


class Gform:
  def __init__(self, in_directory, therm_datafile,pcet=True):
    print therm_datafile
    therm_col = ['mol_tag','E','dHf0','ddH-ddH_ss','H_zpe','S-S_ss']
    self.Gf_df = pd.DataFrame(columns=therm_col)
    self.therm_data = pd.read_csv(therm_datafile, index_col=0) #Contains thermo ref info
    self.indir = in_directory
    #Conversion factors and constants
    self.T = 298.15 #The assumed temperature
    self.Ht2kC = 627.509469 #The conversion factor for Hartree to kCal/mol
    self.vib_scale = 0.9512 #An empirical parameter for fitting DFT freq to experiment
    self.kCal2kJ = 4.184 
    self.Htrpv = 2.372
    #SSmPa is a dictionary with molecules per atom in standard state
    self.SSmPa = {'H':0.5,'C':1.0,'N':0.5,'O':0.5,'F':0.5,'Ni':1.0,'S':1.0}
    self.pcet = pcet

  def dG(self, infile, fq_file=None, pcm_file=None):
    print "Performing PCET Calc: " + infile
    qdat = Qdata.Qdata(infile)
    if fq_file == None: qdat.loadFq(infile)
    else: qdat.loadFq(fq_file)
    print "frequencies?: ",infile,fq_file
    if qdat.H == None: raise ValueError("Couldn't find necessary thermo info")
    if pcm_file == None: qdat.loadPCM(infile); print "No PCM file found, PCM infile?";
    else: qdat.loadPCM(pcm_file); print "Found independent PCM file";
    dGfaq = self.calcG(qdat) #overwrites qdat.E with E including solv effects
    return dGfaq

  def calcG(self, qdat):
    Eatoms, dHfatoms, ddH_ss, S_ss = self.thermSS(qdat)
    ddH = qdat.Hvib*self.vib_scale - qdat.Hzpe*self.vib_scale + self.Htrpv
    dHAtomization = -1.0*qdat.E*self.Ht2kC + Eatoms*self.Ht2kC - qdat.Hzpe*self.vib_scale
    dHf0 = dHfatoms - dHAtomization
    dHf298 = dHf0 + ddH - ddH_ss
    dGfg = dHf298 - self.T*(qdat.S - S_ss)
    print 'E: ' + str(qdat.E)
    print 'dHf0: ' + str(dHf0)
    print 'ddH - ddH_ss: ' + str(ddH - ddH_ss)
    print 'H_zpe: ' + str(qdat.Hzpe)
    print 'S-S_ss: ' + str(self.T*(qdat.S - S_ss))
    dGfaq = dGfg + 1.89 #+ qdat.Gpcm #qdat.loadPCM overwrites qdat.E w/ E including solv effects
    return [name,qdat.E,dHf0,ddH-ddH_ss,qdat.Hzpe,qdat.S-S_ss,dGfaq]
#    return dGfaq
#    return dHf298

  def dGChem(self, infile, pcm_file):
    print 'and you thought you had won!'
    qdat = Qdata.Qdata(infile)
    if qdat.H == None:
      print infile + " is not a frequency calculation."
    qdat.loadPCM(pcm_file)
    Eatoms, dHfatoms, ddH_ss, S_ss = self.thermSS(qdat)
    H298 = qdat.E*self.Ht2kC + qdat.Hvib*self.vib_scale + self.Htrpv
    G = H298 - self.T*qdat.S + qdat.Gpcm + 1.89
    return G

  def batchDG(self):
    for filename in os.listdir(self.indir):
     if filename.split('.')[-1] == 'out':
      fspl = filename.split('_')
      job = fspl[-2]
      if job == 'sfq':
        dG_formation = self.dG(self.indir+"/"+filename,pcm_file=(self.indir+"/"+filename))
      elif job == 'opfqSo' or job == 'sopfqSo':
        chem_id = self.chemID(filename)
        print "getting dG_formation"
        dG_formation = self.dG(self.indir+'/'+filename)
      elif ('freq' in job or 'fq' in job):
        chem_id = self.chemID(filename)
        for pcm_file in os.listdir(self.indir):
          if ('spSo' in pcm_file or 'Sm8' in pcm_file or 'Sm12' in pcm_file or 'sop' in pcm_file):
            pcm_id = self.chemID(pcm_file)
            if chem_id == pcm_id:
              if self.pcet=='True':
                dG_formation = self.dG(self.indir+"/"+filename,pcm_file=(self.indir+"/"+pcm_file))
              else:
                dG_formation = self.dGChem(self.indir+"/"+filename,self.indir+"/"+pcm_file)
      else: continue
      print filename + ": " + str(dG_formation)
      self.Gf_df.loc[len(self.Gf_df)] = 
      print " "

#freq calcs are expensive, this uses data from freq calcs on smaller graphene sheets
  def calledDG(self,fq_dir): #fq_list is a set of base freq files to read
#This part reads the fq_list and pulls thermo data
    therm_subst_dict = {} #dictionary for therm data on different substituents
    for filename in os.listdir(fq_dir):
      true_file = filename.split('/')[-1]
      fspl = true_file.split('_')
      bonded_sub = self.nameChemSub(fspl[0])
      if bonded_sub == None: raise ValueError('Ref substituent unreadable')
      qdat = Qdata.Qdata(fq_dir+'/'+filename)
      qdat.loadFq(fq_dir+'/'+filename)
#      therm_subst_dict[bonded_sub] = [qdat.H,qdat.Hvib,qdat.Hzpe,qdat.S]
      therm_subst_dict[bonded_sub] = qdat
      print therm_subst_dict

#This part reads the real molecule solvent files and calcs dG with corresponding thermo data
    for filename in os.listdir(self.indir):
      if filename.split('.')[-1] == 'out':
        fspl = filename.split('_')
        job = fspl[-2]
        if job == 'opSo' or job == 'spSo':
          print fspl[0]
          bonded_sub = self.nameChemSub(fspl[0])
          if bonded_sub == None:
            print "Couldn't figure out what was bound on: "+filename
            continue
#          therm_subst = therm_subst_dict[bonded_sub]
#          print filename+': '+bonded_sub
#          print therm_subst
          therm_subst = Qdata.Qdata(self.indir+'/'+filename)
          therm_subst.loadPCM(self.indir+'/'+filename)
          therm_subst.readQdat(therm_subst_dict[bonded_sub],'fq')
          dGfaq = self.calcG(therm_subst)
          print filename+": "+str(dGfaq)

#figures out, based on the name, what molecule is bound to the ORR catalyst
  def nameChemSub(self, molecule_name):
    if 'OHs5O2s6' in molecule_name: sub = 'OHs5O2s6'
    if 'OHs5OHs6' in molecule_name: sub = 'OHs5O2s6'
    if 'OHs5O2Hs6' in molecule_name: sub = 'OHs5O2s6'
    if 'OHs5Os6' in molecule_name: sub = 'OHs5O2s6'
    if 'O2H' in molecule_name: sub = 'O2H'
    elif 'O2' in molecule_name: sub = 'O2'
    elif 'OH' in molecule_name: sub = 'OH'
    elif 'O' in molecule_name: sub = 'O'
    else:
        print "happening?"
        if molecule_name[:3] == 'pyc': 
            sub = 'bare' #may revisit this if unreliable
            print "did this even happen?"
#            if len(molecule_name) <= 4: sub = 'bare' #indicates nothing is bound
#            else: sub = None  #This indicates something is wrong
        elif molecule_name[:4] == 'fpda': sub = 'bare' 
#            if len(molecule_name) <= 5: sub = 'bare'
#            else: sub = None
        else: sub = None
    return sub

  def chemID(self, qoutname):
    splile = qoutname.split('.')[0]
    splile = splile.split('_')
    return splile[0] + splile[-1]

  def thermSS(self, Qdata):
    dHfatoms = 0.0
    S_ss = 0.0
    ddH_ss = 0.0
    Eatoms = 0.0
    atom_count = Counter(Qdata.atoms)
    for key_atom in atom_count:
      dHfatoms += atom_count[key_atom] * self.therm_data['H0'][key_atom]/self.kCal2kJ
      ddH_ss += atom_count[key_atom] * self.therm_data['dHss'][key_atom]*\
                self.SSmPa[key_atom]/self.kCal2kJ
      S_ss += atom_count[key_atom] * self.therm_data['S'][key_atom]*self.SSmPa[key_atom]/self.kCal2kJ
      Eatoms += atom_count[key_atom] * self.therm_data['E'][key_atom]
    return Eatoms, dHfatoms, ddH_ss ,S_ss

if __name__ == "__main__":
  in_dir = sys.argv[1]
  pcet = sys.argv[2]
  tdat = "/home/nricke/PyMod/tdat.csv"
  if len(sys.argv) == 3:
    gform = Gform(in_dir,tdat, pcet)
    gform.batchDG()
  elif len(sys.argv) > 3:
    fq_dir = sys.argv[3]
    gform = Gform(in_dir,tdat, pcet)
    gform.calledDG(fq_dir)




