import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import os
import glob, os
#class QchemData


#Generalized set of functions for working with Q-Chem output files

charge_tran = {'a1':-1,'a0':0,'c1':1,'c2':2}
revd=dict([reversed(i) for i in charge_tran.items()])
charge_tran.update(revd)

def sectionGrab(key_start,key_end,infile,start_cut=0,end_cut=0):
  y_read = 0
  grab_list = []
  with open(infile,'r') as f:
    for line in f:
      if key_start in line:
        #[f.next() for i in range(start_cut)] #this is a possible option
        y_read = 1
      if y_read == 1:
        if key_end in line:
          y_read = 0
          break
        else: 
          grab_list.append(line.strip())

  if end_cut == 0:
    return grab_list[start_cut:]
  else:
    return grab_list[start_cut:(-1*end_cut)]

def chargeGrab(in_file,charge=0):
 if charge=="mulliken":
  grab_start = "Ground-State Mulliken Net Atomic Charges"
 elif charge=="chelpg" :
  grab_start = "Ground-State ChElPG Net Atomic Charges"
 else: raise NotImplementedError

 grab_end = "Sum of atomic charges"
 charge_raw = sectionGrab(grab_start,grab_end,in_file,4,1)
 charge = np.zeros(len(charge_raw))
 atom = []
 i = 0
 for line in charge_raw:
  spline = line.split()
  charge[i] = float(spline[2])
  atom.append(spline[1])
  i += 1
 return charge, atom

def coordGrab(in_file):
 grab_start = 'Standard Nuclear Orientation (Angstroms)'
 grab_end = 'Nuclear Repulsion Energy'
 raw_coord = sectionGrab(grab_start,grab_end,in_file,3,1)
 coord = np.zeros((len(raw_coord),3))
 i = 0
 for line in raw_coord:
  spline = line.split()
  coord[i,:] = [float(j) for j in spline[2:]]
  i += 1
 return coord

def optCoordGrab(in_file,out_file=None):
  key_start = "Standard Nuclear Orientation (Angstroms)"
  key_end = "--------------"
  if out_file != None: f_out = open(out_file,'w')

  y_read = False
  with open(in_file,'r') as f:
    for line in f:
      if key_start in line:
        y_read = True
        grab_list = []                    #dump current grab_list
      elif key_end in line:
        y_read = False
        if out_file != None:
          for item in grab_list:
            f_out.write('%s\n' % item)
      elif y_read == True:
        grab_list.append(line.strip())
  raw_coord = grab_list[2:]
  xyz = []
  for item in raw_coord:
    xyz.append(' '.join(item.split()[1:]))
  return xyz

#under development
def freqParse(in_file,in_dir):
  for line in file:
    if "Total Enthalpy" in line:
      pass

def renameFiles(dir, pattern, ext_pattern):
    print dir
    for pathAndFilename in glob.iglob(os.path.join(dir, pattern)):
        print 'Its alive!'
        title, ext = os.path.splitext(os.path.basename(pathAndFilename))
        os.rename(pathAndFilename, 
                  os.path.join(dir, title + ext_pattern))

def SP(fname,coord,charge,mult,method,basis):
 in_file = open(fname,'w')
 in_file.write('molecule\n')
 in_file.write(charge + ' ' + mult)
 np.coord.tofile(in_file)
 in_file.write('$end\n\n $rem')
 in_file.write('method' + ' ' + method + '\n')
 in_file.write('basis' + ' ' + basis + '\n')
 in_file.write('lowdin_population true \n print_orbitals true \n molden_format true \n')
 in_file.write('$end\n')
 in_file.close()

def opt2Sp(in_dir):
 for filename in os.listdir(in_dir):
  functional = 'b3lyp'
  basis = '6-31+g*'
  coord = coordGrab(filename,in_dir)
  sp_fname = filename.split('_')
  ch_name = sp_fname[-1][0:2]
  ch_charge = charge_tran[ch_name]
  an_charge = ch_charge - 1
  an_name = charge_tran[an_charge]
  jobtype = 'sp'
  job1 = sp_fname[0] + '_' + jobtype + '_' + ch_name + '.in'
  job2 = sp_fname[0] + '_' + jobtype + '_' + an_name + '.in'
  SP(job1,coord,ch_charge,1,functional,basis)
  SP(job2,coord,an_charge,2,functional,basis)

#opt2Sp('042715/q_out')

















