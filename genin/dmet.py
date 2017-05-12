import numpy
import logging
kennyloggins = logging.getLogger('dmet')

import schmidt
import xform
import hf
import nr
import partition
import config
import al

from bootstrap import Bootstrap
from dmetenergy import DMETEnergy




class OrbitalOptimized(DMETEnergy):
    def __init__(s, lattice, part, solver):
        s.lat = lattice
        s.solver = solver

        s.hfull,s.U = s.lat.get_operators()

        s.part = part
        s.sites = s.part.lattice_indices()
        s.notsites = list(set(range(len(s.hfull))) - set(s.sites))
        s.ni = len(s.sites)

    def _dmet_iter(s, pot):
        hpot = numpy.zeros(s.hfull.shape)
        for i,item in enumerate(pot):
            hpot[i,i] = item
        #kennyloggins.debug("hpot:\n%s" % str(hpot))

        #w,v = numpy.linalg.eigh(s.hfull + hpot) # apparently the OEP problem 
        #                                    # maps onto the same space
        #                                    # as the HF problem as long
        #                                    # as the potentials and interactions
        #                                    # are on-site
        #psi = v[:,:s.lat.m]
        fock = hf.HF(s.hfull + hpot,s.U,s.lat.m)
        psi = fock.get_Cocc()
        
        T = schmidt.schmidt(psi, s.sites) # fast schmidt?
        himp = xform.one(T, s.hfull)
        Vimp = xform.twoU(T, s.U)

        s.h = himp
        s.V = Vimp

        s.solver.solve(himp, Vimp, s.ni*2)

        target_pop = s.lat.filling()*2 # XXX: factor of two?
        return s.solver.energy, numpy.diag(s.solver.P)[:s.ni]-target_pop

    def optimize(s):
        guess = numpy.zeros(len(s.hfull)/2 + len(s.hfull)%2 - s.ni/2)
        def wrap(pot):
            longpot = numpy.zeros(len(s.hfull))
            for i in range(len(pot)):
                longpot[i+s.ni] = pot[i]
                longpot[-(i+1)] = pot[i]
            for i in s.sites:
                longpot[i] = 0
            print pot
            return s._dmet_iter(longpot)

        optpot = al.al(wrap, guess)
        s._dmet_iter(optpot) # Just to make sure
        s._calc_energy()
        kennyloggins.debug("Final density matrix:\n"+str(s.solver.P))

    def opt_cons(s):
        guess = numpy.zeros(s.hfull.shape[0]-len(s.sites))
        optpot = nr.nr(lambda x: s._dmet_iter(x)[1], guess)
        s._dmet_iter(optpot) # Just to make sure
        s._calc_energy()
        kennyloggins.debug("Final density matrix:\n"+str(s.solver.P))



class Increments:
    def __init__(s, lattice, solver, order=2):
        s.lat = lattice
        s.solver = solver

        assert(lattice.ny == 1)

        # Compute HF for system - we only need to do this once, but we'd have to refactor bootstrap
        #h,U = s.lat.get_operators()
        #fock = hf.HF(h,U,s.lat.m)
        #s.psi = fock.get_Cocc()

        s.e1 = None
        s.e2 = None
        s.e = None


    def optimize(s):
        kennyloggins.info("Beginning method-of-increments bootstrap calculation")
        n = s.lat.nsites()


        # one site part
        part1s = partition.Partition(s.lat,[[0,0]],[[0]],[0],[])
        kennyloggins.info("Running 1 site calculation")
        bs1s = Bootstrap(s.lat, part1s, s.solver)
        bs1s.optimize()
        #e1s = bs1s.energy()
        e1s = bs1s.eci

        #two site part
        e2s = numpy.zeros(n)
        e2s[0] = 2*e1s
        for j in range(1,n):
            kennyloggins.info("Running 2 site calculation %d/%d" % (j,n-1))
            part2s = partition.Partition(s.lat,[[0,0],[j,0]],[[0,1]],[0,1],[])
            bs2s = Bootstrap(s.lat, part2s, s.solver)
            bs2s.optimize()
            #e2s[j] = bs2s.energy()*2 # factor of 2 because the reported energy is persite
            e2s[j] = bs2s.eci
        kennyloggins.info("singles energy: %f" % e1s)
        kennyloggins.info("pair energies:  %s" % str(e2s))
        numpy.set_printoptions(precision=6)
        kennyloggins.info("pair increment energy:  %s" % str(e2s-2*e1s))

        efull = 0
        for i in range(n):
            efull += e1s 
        for i in range(n):
            for j in range(1,n):
                efull += .5*(e2s[j] - e1s - e1s)
        efull /= n
                

        s.e = (-n + 2)*e1s + .5*e2s[1:].sum()
        kennyloggins.debug("reduced energy: %f" % s.e)
        kennyloggins.debug("full energy:    %f" % efull)
        s.e = efull

    def energy(s):
        return s.e
    def center_energy(s):
        return 0




