#!/usr/bin/env python
import os
import subprocess


"""
Script for running hringDmrg.py. The pyscf script needs to run in the same directory where it produces its output at the moment, ie FCIDUMP is generated in the cwd and then the output need to be there as well. I suppose this means I need to either change the path or find some way to change where it looks for the FCIDUMP file.
"""

#st_rings = [str(i) for i in range(8,23,2)]
st_rings = [str(i) for i in range(8,13,2)]

time = '2-00'

hring_code = '/home/nricke/PyMod/pyscfScripts/hringDmrg.py'
init_dir = os.getcwd()

#run code in a directory that we want to populate with subdirectories of jobs
for hring in st_rings:
    job_name = 'h' + hring
    jobdir = "%s/%s" % (init_dir,job_name)
    outfile = hring+'.out'
    errfile = hring+'.err'
    #make subdirectory
    if not os.path.exists(jobdir):
        os.makedirs(jobdir)
    os.chdir(jobdir) #change to subdirectory
    job = ['sqthis','-t',time,'-J',job_name,'--mem','2000','-o',outfile,'-e',errfile,hring_code,hring]
    subprocess.call(job)
    os.chdir(init_dir)


