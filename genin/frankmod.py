#! /usr/bin/env python
import argparse
import numpy
import logging
import sys

import fragment
import bootstrap
import nr
import dmetenergy

kennyloggins = logging.getLogger('frankentonian')

class Frankentonian(dmetenergy.DMETEnergy):
    def __init__(s, lat):
        s.lat = lat

        #Build all pairs fragments
        #TODO: Update for non-translationally awesome stuff
        s.one_frag = fragment.Fragment("0")
        s.one_frag.add_site("A", 1)
        s.one_frag.center=["A"]
        s.one_frag.index()

        s.two_frags = []
#        for j in range(1, s.lat.n/2+1):
#            frag = fragment.Fragment("0%d" % j)
#            frag.add_site("A",0)
#            frag.add_site("B",j)
#            frag.index()
#            s.two_frags.append(frag)
#            frag = fragment.Fragment("0%d" % -j)
#            frag.add_site("A",0)
#            frag.add_site("B",-j)
#            frag.index()
#            s.two_frags.append(frag)
#            if j == 1: break

#        frag = fragment.Fragment("0%d" % 1)
#        frag.add_site("A",0)
#        frag.add_site("B",1)
#        frag.add_site("C",2)
#        frag.center=["A"]
#        frag.index()
#        s.two_frags.append(frag)
#
        frag = fragment.Fragment("0%d" % 1)
        frag.add_site("A",0)
        frag.add_site("B",1)
        frag.index()
        s.two_frags.append(frag)

#        frag = fragment.Fragment("0%d" % 3)
#        frag.add_site("A",1)
#        frag.add_site("B",2)
#        frag.index()
#        s.two_frags.append(frag)

        print s.two_frags

        #Build constraints
        constraints = {}
        constraints['filling'] = []
        constraints['1e'] = []
        constraints['2e'] = []

        checks = {}
        checks['filling'] = []
        checks['1e'] = []
        checks['2e'] = []

        # Filling constraints and checks
        constraints['filling'].append([
            {'fragment': s.one_frag.name,
             'site'    : s.one_frag.sites[0]['name']}])
        checks['filling'].append(
            {'fragment': s.one_frag.name,
             'site'    : s.one_frag.sites[0]['name']})

        for frag in s.two_frags:
            constraints['filling'].append([
                {'fragment': frag.name,
                 'site'    : frag.sites[0]['name']},
                {'fragment': frag.name,
                 'site'    : frag.sites[1]['name']}
                ])

            checks['filling'].append(
                {'fragment': frag.name,
                 'site'    : frag.sites[0]['name']})

            checks['filling'].append(
                {'fragment': frag.name,
                 'site'    : frag.sites[1]['name']})

        # Two electron pop stuff
        good = {'fragment': s.one_frag.name,
                'site'    : [s.one_frag.sites[0]['name']]*4}
        for frag in s.two_frags:
            bad = []
            for site in frag.sites:
                bad.append( 
                        {'fragment': frag.name,
                        'site'    : [site['name']]*4})
            
            constraint = {'good': good, 'bad': bad}
            constraints['2e'].append(constraint)

            for i in range(len(bad)):
                check = {'good': good,
                         'bad' : bad[i]}
                checks['2e'].append(check)


        s.part = partition.Partition([s.one_frag]+s.two_frags, constraints, checks)
        #Constuct bootstrap object
        s.bs = bootstrap.Bootstrap(s.lat, s.part)

    def _frank_iter(s, pot):
        # Build bootstrap potentials
        s.bs._make_pots(pot)
        # Solve each two site fragment (to get Hamiltonian increment)
        # Also, solve the one site fragment as a hack to get its Hamiltonian
        for frag in [s.one_frag] + s.two_frags:
#        for frag in s.two_frags:
            frag.solver.solve(frag.h+frag.hpot,frag.V+frag.Vpot,hnat=frag.h,Vnat=frag.V,nelec=frag.ni)
#            print "two",frag.solver.get_P2()[0,0,0,0]
#            print frag.solver.get_P()[0,0]
            print "Fragment:", frag.lattice_indices()
            print frag.solver.H_red

        # Build the method of increments hamiltonian
#        HA = s.one_frag.solver.H_red # FCI code should reduce 1-site hamiltonian to itself
#        print HA
        HA_dressed = s.H0.copy() #HA.copy()
        for frag in s.two_frags:
            # using more verbose form so we can check size of increment
#            print "Fragment:", frag.lattice_indices()
#            print "H_reduced"
#            print frag.solver.H_red
            increment = frag.solver.H_red - s.H0
            
            # if good approx., increment should be much smaller than HA
            if numpy.linalg.norm(increment) > numpy.linalg.norm(s.H0):
                kennyloggins.warning("Increment is larger than H0!")
            HA_dressed += increment
#            print increment

        try:
            w,v = numpy.linalg.eigh(HA_dressed)
        except numpy.linalg.linalg.LinAlgError as err:
            print "Linalg error while diagonalizing HA_dressed in Frank Iter; print HA_dressed"
            print HA_dressed
            print "HA:"
            print HA
            sys.exit(-1)
        psi = v[:,0] # ground state
        psi[2] = -psi[2] # I think we have a sign problem
        psi[1] = -psi[1] # relative to Troy
        X = psi.reshape((2,2))

        P, G = s.one_frag.solver.calcPG(2,1,X)

        #pretty hacky
        s.one_frag.solver.P[:,:,0,0] = P[:,:]
        s.one_frag.solver.P2 = G

#        print s.one_frag.solver.get_P()[0,0]
        print "one",s.one_frag.solver.get_P2()[0,0,0,0]

        obj = s.bs._calc_obj()
        return obj


    def optimize(s):
        guess = numpy.zeros(s.part.optlen())
        s.one_frag.solver.solve(s.one_frag.h,s.one_frag.V,nelec=s.one_frag.ni)
        s.H0 = s.one_frag.solver.H_red
        print "H0:"
        print s.H0
        obj = s._frank_iter(guess) # Just to make sure
#        sys.exit(-1)

        optpot = nr.nr(s._frank_iter, guess)
        s._frank_iter(optpot) # Just to make sure
        s.bs._check()

        s.bs._calc_energy(s.one_frag)
        s.e = s.one_frag.e
        s.ec = s.one_frag.ec
        s.e1 = s.one_frag.e1
        s.e2 = s.one_frag.e2

        print "One-site fragment:"
#        print s.one_frag.solver.get_P()
        print s.one_frag.solver.get_P2()[0,0,0,0]

        for frag in s.two_frags:
            print "Fragment:", frag.lattice_indices()
            print "P:",frag.solver.get_P()[0,1]
#            print frag.solver.get_P()
            print "P2:"
            diag = []
            P2 = frag.solver.get_P2()
#            for i in range(P2.shape[0]):
#                diag.append(P2[i,i,i,i])
#            print diag
            print P2[0,0,0,0]
            print

#        print "vA eigenvector: ",s.vA0
#        print "dressed eigenvector: ",s.vA0d

        return optpot







if __name__ == "__main__":
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    root.addHandler(ch)

    import lattice
    import solvers
    import config
    import partition

    numpy.set_printoptions(precision = 3, suppress = True, linewidth=20000)

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', help='Number of electron pairs', type=int, required=True)
    parser.add_argument('-l', '--lattice',  help='Lattice file', type=str, required=True)
    parser.add_argument('--core', help='Use core Hamiltonian correction', type=bool, default=True)
    parser.add_argument('--solver', help="Impurity solver (GK|Troy)", type=str, default="Troy")
    parser.add_argument('--debug', help="Debug output", action='store_true')
    parser.add_argument('--no-two', help='Disable two electron checks and constraints', action='store_true')
    args = parser.parse_args()

    if args.debug:
        ch.setLevel(logging.DEBUG)

    kennyloggins.debug("Input:\n%s" % args)

    # Update globals from input
    config.bootstrap['core'] = args.core
    config.bootstrap['no_two'] = args.no_two
    config.solvers['type'] = args.solver

    l = lattice.Lattice.fromfile(args.lattice, args.m)

    kennyloggins.info("Performing frankentonian optimization")

    # Call frankentonian
    emb = Frankentonian(l)
    emb.optimize()

    kennyloggins.info("Final energy: %f" % emb.energy())
    kennyloggins.info("Center energy: %f" % emb.center_energy())

    kennyloggins.info("One electron energy: %f" % emb.e1)
    kennyloggins.info("Two electron energy: %f" % emb.e2)

    kennyloggins.info("Done")
    logging.shutdown()
