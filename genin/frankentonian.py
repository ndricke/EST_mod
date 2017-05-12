import argparse
import numpy
import logging
import sys

import hf
import xform
import schmidt

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
        for j in range(1, s.lat.n):
#        for j in range(1,2):
            frag = fragment.Fragment("0%d" % j)
            frag.add_site("A",0)
            frag.add_site("B",j)
#            frag.add_site("C",j+1)
            frag.index()
            s.two_frags.append(frag)

#            frag = fragment.Fragment("0%d" % (-j))
#            frag.add_site("A",0)
#            frag.add_site("B",-j)
#            frag.index()
#            s.two_frags.append(frag)
#        print s.two_frags

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

        for frag in s.two_frags:
            site_A = frag.name_to_lat_idx['A']
            site_B = frag.name_to_lat_idx['B']
            frag.Ihfull,frag.Ivfull = s.generateI(site_A,site_B,s.lat.h,s.lat.V)

            hc1 = xform.core(s.one_frag.Tc,frag.Ivfull)
            Ih1 = xform.one(s.one_frag.T,frag.Ihfull+hc1)
            Iv1 = xform.two(s.one_frag.T,frag.Ivfull)

            hc2 = xform.core(frag.Tc,frag.Ivfull)
            frag.Ih2 = xform.one(frag.T,frag.Ihfull+hc2)
            frag.Iv2 = xform.two(frag.T,frag.Ivfull)

            X1 = numpy.zeros((2,2)); X1[0,0] = 1.0
            X2 = numpy.zeros((2,2)); X2[1,0] = 1.0
            X3 = numpy.zeros((2,2)); X3[0,1] = 1.0
            X4 = numpy.zeros((2,2)); X4[1,1] = 1.0
            frag.I0 = s.one_frag.solver.calcHred(1,Ih1,Iv1,X1,X2,X3,X4)

#site_1 is the single site we project fully down to, site_2 is the interim auxiliary site
    def generateI(s,site_1,site_2,h,V):
            n = h.shape[0]
            Ih = numpy.zeros((n,n))
            Iv = numpy.zeros((n,n,n,n))
#I put the 0.5 terms first so they are overwritten by the A-B interactions            
            Ih[:,site_2] = 0.5*h[:,site_2]
            Ih[site_2,:] = 0.5*h[site_2,:]
            Ih[site_1,site_2] = h[site_1,site_2]
            Ih[site_2,site_1] = h[site_2,site_1]

            Iv[site_2,:,site_2,:] = 0.5*V[site_2,:,site_2,:]
            Iv[:,site_2,:,site_2] = 0.5*V[:,site_2,:,site_2]
            Iv[site_1,site_1,site_1,site_1] = V[site_1,site_1,site_1,site_1]
            Iv[site_2,site_2,site_2,site_2] = V[site_2,site_2,site_2,site_2]

#            print Ih
#            for i in range(n):
#                print Iv[i,i,i,i]

            return Ih,Iv

    def _frank_iter(s, pot):
        # Build bootstrap potentials
        s.bs._make_pots(pot)

        # Solve each two site fragment (to get Hamiltonian increment)
        s.one_frag.solver.solve(s.one_frag.h+s.one_frag.hpot,s.one_frag.V+s.one_frag.Vpot,hnat=s.one_frag.h,Vnat=s.one_frag.V,nelec=s.one_frag.ni)

        HA = s.one_frag.solver.H_red # FCI code should reduce 1-site hamiltonian to itself
        HA_dressed = HA.copy()
        print HA_dressed
        for frag in s.two_frags:
            frag.solver.solve(frag.h+frag.hpot,frag.V+frag.Vpot,hnat=frag.Ih2,Vnat=frag.Iv2,nelec=frag.ni)
            increment = frag.solver.H_red - frag.I0
            #if this is a good approx, the increment should be smaller than HA
            if numpy.linalg.norm(increment) > numpy.linalg.norm(HA):
                kennyloggins.warning("Increment is larger than HA!")
            HA_dressed += increment
            print "H_red:"
            print frag.solver.H_red
            print "Increment:"
            print increment

        print HA_dressed

        # Diagonalize dressed one-site hamiltonian and build density matrices
        w,v = numpy.linalg.eigh(HA_dressed)
        psi = v[:,0] # ground state
        psi[2] = -psi[2] # I think we have a sign problem
        psi[1] = -psi[1] # relative to Troy
        X = psi.reshape((2,2))

        P, G = s.one_frag.solver.calcPG(2,1,X)

        #pretty hacky
        s.one_frag.solver.P[:,:,0,0] = P[:,:]
        s.one_frag.solver.P2 = G

        obj = s.bs._calc_obj()
        return obj


    def optimize(s):
        guess = numpy.zeros(s.part.optlen())

        obj = s._frank_iter(guess) # Just to make sure
        sys.exit(-1)

        optpot = nr.nr(s._frank_iter, guess)
        s._frank_iter(optpot) # Just to make sure
        s.bs._check()

        s.bs._calc_energy(s.one_frag)
        s.e = s.one_frag.e
        s.ec = s.one_frag.ec
        s.e1 = s.one_frag.e1
        s.e2 = s.one_frag.e2

        print "One-site fragment:"
        print s.one_frag.solver.get_P()
        print s.one_frag.solver.get_P2()[0,0,0,0]

        for frag in s.two_frags:
            print "Fragment:", frag.lattice_indices()
            print "P:"
            print frag.solver.get_P()
            print "P2:"
            diag = []
            P2 = frag.solver.get_P2()
            for i in range(P2.shape[0]):
                diag.append(P2[i,i,i,i])
            print diag
            print

        return optpot

    def hubbImpE(s,h,V,P1,P2):
        e1 = 0
        e2 = 0
        n = h.shape[0]
        ni = n//2

        for f in range(ni):
            for j in range(2*ni):
                e1 += h[f,j] * P1[j,f]

        for f in range(ni):
            e2 += V[f,f,f,f]*P2[f,f,f,f]

        e1 /= ni
        e2 /= ni
        e = e1 + e2
        return e,e1,e2

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

    kennyloggins.debug("One electron energy: %f" % emb.e1)
    kennyloggins.debug("Two electron energy: %f" % emb.e2)

    kennyloggins.info("Done")
    logging.shutdown()
