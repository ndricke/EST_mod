#!/usr/bin/env python
import os
import subprocess

header = "#!/bin/sh"

sbatch_list = """
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=2000
#SBATCH -x n03
#SBATCH -x h05
"""
#SBATCH -o h8.out
#SBATCH -e h8.err

export_lines = """
export PYTHONPATH=${PYTHONPATH}:/home/nricke/.local/lib/python2.7/site-packages/pip/_vendor
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/mmavros/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/pdesilva/software/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64_lin
"""

st_rings = [str(i) for i in range(8,31,2)]
#st_rings = [str(i) for i in range(8,9,2)]

time = '2-00'
sbatch_file = "submit.sh"

hring_code = '/home/nricke/PyMod/pyscfScripts/hringDmrg.py'
init_dir = os.getcwd()

#run code in a directory that we want to populate with subdirectories of jobs
for n in st_rings:
    job_name = 'h' + n
    jobdir = "%s/%s" % (init_dir,job_name)
#    outfile = n+'.out'
#    errfile = n+'.err'
    #make subdirectory
    if not os.path.exists(jobdir):
        os.makedirs(jobdir)
    os.chdir(jobdir) #change to subdirectory
    with open(sbatch_file,'w') as f:
        f.write(header)
        f.write(sbatch_list)
#        f.write(sbatch_output)
        f.write(export_lines)
        f.write("out_path=./\n")
#        f.write("cd ${out_path}\n")
        f.write("~/PyMod/pyscfScripts/hringDmrg.py %s > ${out_path}/job.log | tee ${out_path}/err.log\n" % n)

#    job = ['sqthis','-t',time,'-J',job_name,'--mem','2000','-o',outfile,'-e',errfile,hring_code,hring]
    job = ['sbatch','-J',job_name,sbatch_file]
    subprocess.call(job)
    os.chdir(init_dir)


