#!/bin/sh

#SBATCH -t 480:00:00
#SBATCH -o /home/hzye2011/projects/carbon_dimer/first_trial/cc-pvqz/atom/1DS/sbatch.out
#SBATCH -e /home/hzye2011/projects/carbon_dimer/first_trial/cc-pvqz/atom/1DS/sbatch.err
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem=20000
#SBATCH -x n03
#SBATCH -x h05

# export lib path for user specific lib
export PYTHONPATH=${PYTHONPATH}:/home/hzye2011/.local/lib/python2.7/site-packages/pip/_vendor
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/mmavros/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/pdesilva/software/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64_lin

out_path=/home/hzye2011/projects/carbon_dimer/first_trial/cc-pvqz/atom/1DS

cd ${out_path}
python c.py 2> ${out_path}/err.log | tee ${out_path}/job.log
#rm -rf /scratch/hzye2011/dmrg/atom_1DS-qz
