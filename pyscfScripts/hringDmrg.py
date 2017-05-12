#!/usr/bin/env python
import numpy as np
import sys
import os
import hringGen
import argparse
from pyscf import gto, scf, dmrgscf, mcscf


def runDmrg(n,m,in_atom,in_basis='sto-3g'):
    mol = gto.M(atom=in_atom, basis=in_basis)
    mf = scf.RHF(mol)
    print mf.kernel()

    #mc = dmrgscf.DMRGSCF(mf, 8, 8)
    #print mc.kernel()

    mc = mcscf.CASCI(mf, n, m) #should be about the same?
    mc.fcisolver = dmrgscf.DMRGCI(mol)
    mc.fcisolver.spin = 0

    #   setting for output
    out_path= os.getcwd()
    mc.fcisolver.runtimeDir = out_path
    mc.fcisolver.scratchDirectory = "/scratch/nricke/"+out_path.split('/')[-1] #name scratch after job
    mc.fcisolver.memory = 2

    #   setting for run
    mc.fcisolver.mpiprefix = "mpirun -np 1"

    #   schedule
    mc.fcisolver.scheduleSweeps = [0, 4, 8, 12, 16, 20, 24, 30]
    mc.fcisolver.scheduleMaxMs  = [200, 400, 800, 1200, 2000, 2000, 2000, 2000]
    mc.fcisolver.scheduleTols   = [0.0001, 0.0001, 0.0001, 0.0001, 1e-5, 1e-7, 1e-7, 1e-7]
    mc.fcisolver.scheduleNoises = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0, 0.0]
    mc.fcisolver.twodot_to_onedot = 34
    mc.fcisolver.maxIter = 50

    mc.kernel()
    mc.analyze()


def runDmrgHring(n,Hbond=0.951):
    poly = hringGen.InscPoly(n,Hbond)
    in_atom = ' '.join(['H '+' '.join(i)+';' for i in poly.geom])
    runDmrg(n,n,in_atom)

if __name__ == "__main__":
#    parser = argparse.ArgumentParser()
#    parser.add_argument('-n', help='Atoms in H ring', type=str, default='20')
#    parser.add_argument('-d', help='Distance b/w atoms in H ring', type=float, default=0.951)
#    args = parser.parse_args()
#    runDmrgHring(args.n,args.d)

    n = int(sys.argv[1])
    d = float(sys.argv[2])
    runDmrgHring(n,d)


