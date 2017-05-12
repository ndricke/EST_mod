import numpy
import scipy.optimize
import logging
kennyloggins = logging.getLogger('boot')
import itertools
import copy
import sys

import schmidt
import xform
import hf
import nr
import partition
import config
import al

from dmetenergy import DMETEnergy


class Bootstrap(DMETEnergy):
    def __init__(s, lattice, part):
        if not config.bootstrap['core']:
            kennyloggins.warning("Not using core potential!")

        s.lat = lattice
        s.part = part

        #Check that all fragments can be filled
        for frag in s.part.fragments:
            if frag.ni > s.lat.m:
                kennyloggins.error("More sites (%d) in fragment %s than electron pairs in lattice (%d)" % (frag.ni, frag.name, s.lat.m))
                import sys; sys.exit(1)

        # Compute HF for system
        h = s.lat.h.copy()
        V = copy.copy(s.lat.V)
        kennyloggins.info("Beginning HF")
        kennyloggins.info(config.hf)
        s.fock = hf.HF(h,V,s.lat.m,rho_guess=config.hf)

        # Schmidt decompose for each fragment
        s.psi = s.fock.get_Cocc()
        kennyloggins.info("Finished HF")
        kennyloggins.info("HF energy: %f" % s.fock.get_energy()) 
        for fragment in part.fragments:
            sites = fragment.lattice_indices()
            ni = len(sites)

            T = schmidt.schmidt(s.psi,sites)
            Tfull = schmidt.schmidt(s.psi, sites, allv=True)
            Tc = Tfull[:,2*ni:]

            if config.bootstrap['core']:
                # get environment orbitals
                hc = xform.core(Tc,V)
                for i in range(ni): #XXX Why is this step necessary?
                    hc[i,i] *= 0.5
            fragment.h = xform.one(T, h + (hc if config.bootstrap['core'] else 0))
            fragment.V = xform.two(T, V)
            fragment.T = T
            fragment.Tc = Tc
            fragment.Tfull = Tfull

        s.e1 = None
        s.e2 = None
        s.e = None
        print s.part.constraints['filling']

    def _make_pots(s, pot):
        shift = 0
        for frag in s.part.fragments:
            ni = len(frag.sites)
            frag.hpot = numpy.zeros((2*ni, 2*ni))
            frag.Vpot = numpy.zeros((2*ni, 2*ni, 2*ni, 2*ni))

        # Filling
        for con in s.part.constraints['filling']:
            for site in con:
                frag = s._get_fragment_by_name(site['fragment'])
                i = frag.get_site_frag_idx(site['site'])
                frag.hpot[i,i] += pot[shift]
            shift += 1

        # 1e- 
        for con in s.part.constraints['1e']:
            for site in con['bad']:
                frag = s._get_fragment_by_name(site['fragment'])
                i,j = map(frag.get_site_frag_idx, site['site'])
                add_pot = numpy.zeros((frag.ni*2,frag.ni*2))
                add_pot[i,j] = pot[shift]
                add_pot[j,i] = pot[shift]
                frag.hpot += add_pot
            shift += 1

        # 2e- 
        for con in s.part.constraints['2e']:
            for site in con['bad']:
                frag = s._get_fragment_by_name(site['fragment'])
                i,j,k,l = map(frag.get_site_frag_idx, site['site'])
                add_pot = numpy.zeros((frag.ni*2,frag.ni*2,frag.ni*2,frag.ni*2))
                add_pot[i,j,k,l] = pot[shift]
                add_pot[j,i,l,k] = pot[shift]
                add_pot[k,l,i,j] = pot[shift]
                add_pot[l,k,j,i] = pot[shift]
                add_pot[k,j,i,l] = pot[shift]
                add_pot[l,i,j,k] = pot[shift]
                add_pot[i,l,k,j] = pot[shift]
                add_pot[j,k,l,i] = pot[shift]
                frag.Vpot += add_pot
            shift += 1

        # 2e- diff
        if '2e_diff' in s.part.constraints:
            for con in s.part.constraints['2e_diff']:
                for site in con['plus']:
                    frag = s._get_fragment_by_name(site['fragment'])
                    i,j,k,l = map(frag.get_site_frag_idx, site['site'])
                    add_pot = numpy.zeros((frag.ni*2,frag.ni*2,frag.ni*2,frag.ni*2))
                    add_pot[i,j,k,l] = pot[shift]
                    add_pot[j,i,l,k] = pot[shift]
                    add_pot[k,l,i,j] = pot[shift]
                    add_pot[l,k,j,i] = pot[shift]
                    add_pot[k,j,i,l] = pot[shift]
                    add_pot[l,i,j,k] = pot[shift]
                    add_pot[i,l,k,j] = pot[shift]
                    add_pot[j,k,l,i] = pot[shift]
                    frag.Vpot += add_pot
                for site in con['minus']:
                    frag = s._get_fragment_by_name(site['fragment'])
                    i,j,k,l = map(frag.get_site_frag_idx, site['site'])
                    add_pot = numpy.zeros((frag.ni*2,frag.ni*2,frag.ni*2,frag.ni*2))
                    add_pot[i,j,k,l] = pot[shift]
                    add_pot[j,i,l,k] = pot[shift]
                    add_pot[k,l,i,j] = pot[shift]
                    add_pot[l,k,j,i] = pot[shift]
                    add_pot[k,j,i,l] = pot[shift]
                    add_pot[l,i,j,k] = pot[shift]
                    add_pot[i,l,k,j] = pot[shift]
                    add_pot[j,k,l,i] = pot[shift]
                    frag.Vpot -= add_pot
                shift += 1

        # 2e- Energy
        if '2e_energy' in s.part.constraints:
            for con in s.part.constraints['2e_energy']:
                for site in con['bad']:
                    frag = s._get_fragment_by_name(site['fragment'])
                    f = frag.get_site_frag_idx(site['site'])
                    add_pot = numpy.zeros((frag.ni*2,frag.ni*2,frag.ni*2,frag.ni*2))
                    ri = range(frag.ni)

                    # XXX: 
                    # 1) This might not be the right expression when we have exchange
                    #    i.e. it might need to be <fj||kl> instead of <fj|kl>
                    #    (or 2<fj|kl> - <fj|lk>?)
                    # 2) I'm assuming the 2e- integrals have been properly symmeterized 
                    # 3) I'm assuming that the fragment indices are the first ni. This 
                    #    depends on my implementation of Fragment. This should probably
                    #    iterate over the map(frag.get_site_frag_idx, frag.sites)
                    for j,k,l in itertools.product(ri,ri,ri):
                        add_pot[f,j,k,l] = frag.V[f,j,k,l]*pot[shift]
                        add_pot[j,f,l,k] = frag.V[f,j,k,l]*pot[shift]
                        add_pot[k,l,f,j] = frag.V[f,j,k,l]*pot[shift]
                        add_pot[l,k,j,f] = frag.V[f,j,k,l]*pot[shift]
                        add_pot[k,j,f,l] = frag.V[f,j,k,l]*pot[shift]
                        add_pot[l,f,j,k] = frag.V[f,j,k,l]*pot[shift]
                        add_pot[f,l,k,j] = frag.V[f,j,k,l]*pot[shift]
                        add_pot[j,k,l,f] = frag.V[f,j,k,l]*pot[shift]
                    frag.Vpot += add_pot
                shift += 1

        if '2RDM_norm' in s.part.constraints and len(s.part.constraints['2RDM_norm'])>0:
            con = s.part.constraints['2RDM_norm'][0]
            for site in con:
                frag = s._get_fragment_by_name(site['fragment'])
                i,j,k,l = map(frag.get_site_frag_idx, site['site'])
                add_pot = numpy.zeros((frag.ni*2,frag.ni*2,frag.ni*2,frag.ni*2))
                add_pot[i,j,k,l] = pot[shift]
                add_pot[j,i,l,k] = pot[shift]
                add_pot[k,l,i,j] = pot[shift]
                add_pot[l,k,j,i] = pot[shift]
                add_pot[k,j,i,l] = pot[shift]
                add_pot[l,i,j,k] = pot[shift]
                add_pot[i,l,k,j] = pot[shift]
                add_pot[j,k,l,i] = pot[shift]
                frag.Vpot += add_pot
            shift += 1

        assert(shift == len(pot))

    def _calc_obj(s):
        shift = 0
        obj = numpy.zeros(s.part.optlen())

        # Filling
        target = s.lat.filling()*2
        for con in s.part.constraints['filling']:
            site = con[0]
            frag = s._get_fragment_by_name(site['fragment'])
            i = frag.get_site_frag_idx(site['site'])
            obj[shift] = frag.solver.get_P()[i,i] - target
            shift += 1

        # 1e- 
        for con in s.part.constraints['1e']:
            good_site = con['good']
            good_frag = s._get_fragment_by_name(good_site['fragment'])
            i,j = map(good_frag.get_site_frag_idx, good_site['site'])
            good = good_frag.solver.get_P()[i,j]

            bad_site = con['bad'][0]
            bad_frag = s._get_fragment_by_name(bad_site['fragment'])
            i,j = map(bad_frag.get_site_frag_idx, bad_site['site'])
            bad = bad_frag.solver.get_P()[i,j]

            obj[shift] = bad - good
            shift += 1

        # 2e- 
        for con in s.part.constraints['2e']:
            good_site = con['good']
            good_frag = s._get_fragment_by_name(good_site['fragment'])
            i,j,k,l = map(good_frag.get_site_frag_idx, good_site['site'])
            good = good_frag.solver.get_P2()[i,j,k,l]

            bad_site = con['bad'][0]
            bad_frag = s._get_fragment_by_name(bad_site['fragment'])
            i,j,k,l = map(bad_frag.get_site_frag_idx, bad_site['site'])
            bad = bad_frag.solver.get_P2()[i,j,k,l]

            obj[shift] = bad - good
            shift += 1

        # 2e- diff (similar to 2e)
        if '2e_diff' in s.part.constraints:
            for con in s.part.constraints['2e_diff']:
                plus_site = con['plus'][0]
                plus_frag = s._get_fragment_by_name(plus_site['fragment'])
                i,j,k,l = map(plus_frag.get_site_frag_idx, plus_site['site'])
                plus = plus_frag.solver.get_P2()[i,j,k,l]

                minus_site = con['minus'][0]
                minus_frag = s._get_fragment_by_name(minus_site['fragment'])
                i,j,k,l = map(minus_frag.get_site_frag_idx, minus_site['site'])
                minus = minus_frag.solver.get_P2()[i,j,k,l]

                obj[shift] = plus-minus
                shift += 1

        # 2e- Energy 
        if '2e_energy' in s.part.constraints:
            for con in s.part.constraints['2e_energy']:
                good_site = con['good']
                good_frag = s._get_fragment_by_name(good_site['fragment'])
                good_frag_copy = copy.deepcopy(good_frag)
                good_frag_copy.center = [good_site['site']]
                good_frag_copy.index()

                bad_site = con['bad'][0]
                bad_frag = s._get_fragment_by_name(bad_site['fragment'])
                bad_frag_copy = copy.deepcopy(bad_frag)
                bad_frag_copy.center = [bad_site['site']]
                bad_frag_copy.index()

                # XXX: Hijack energy computation for objective function
                #      It would be best to refactor DMETEnergy into 
                #      utility functions instead of a class, I think.
                dmete = DMETEnergy()
                dmete.lat = s.lat
                dmete._calc_energy(good_frag_copy)
                dmete._calc_energy(bad_frag_copy)

                good_e2c = good_frag_copy.e2c
                bad_e2c  = bad_frag_copy.e2c

                obj[shift] = bad_e2c - good_e2c
                shift += 1

        # 2RDM Normalization constraint
        # XXX: this depends on translational symmetry
        # TODO: implement get_drace_gamma 

        if '2RDM_norm' in s.part.constraints and len(s.part.constraints['2RDM_norm'])>0: 
            gamma = get_drace_gamma(s.fragments)
            gammma_ijij = 0 
            n,m = s.lat.n, s.lat.m
            for i,j in itertools.product(range(n), range(n)):
                gamma_ijij += gamma[i,j,i,j]

            obj[shift] = gamma_ijij - (2*m-1)*m
            shift += 1

        assert(shift == len(obj))
        return obj 

        assert(shift == len(obj))
        return obj

    def _check(s):
        _is_zero = lambda x: abs(x)<config.checks['tol']

        # Filling
        target = s.lat.filling()*2
        for check in s.part.checks['filling']:
            site = check
            frag = s._get_fragment_by_name(site['fragment'])
            i = frag.get_site_frag_idx(site['site'])
            if not _is_zero(frag.solver.get_P()[i,i] - target):
                kennyloggins.warning("Failed filling check: %s" % str(check))

        # 1e- 
        for check in s.part.checks['1e']:
            good_site = check['good']
            good_frag = s._get_fragment_by_name(good_site['fragment'])
            i,j = map(good_frag.get_site_frag_idx, good_site['site'])
            good = good_frag.solver.get_P()[i,j]

            bad_site = check['bad']
            bad_frag = s._get_fragment_by_name(bad_site['fragment'])
            i,j = map(bad_frag.get_site_frag_idx, bad_site['site'])
            bad = bad_frag.solver.get_P()[i,j]

            if not _is_zero(bad - good):
                kennyloggins.warning("Failed 1e- check: %s" % str(check))


        # 2e- 
        for check in s.part.checks['2e']:
            good_site = check['good']
            good_frag = s._get_fragment_by_name(good_site['fragment'])
            i,j,k,l = map(good_frag.get_site_frag_idx, good_site['site'])
            good = good_frag.solver.get_P2()[i,j,k,l]

            bad_site = check['bad']
            bad_frag = s._get_fragment_by_name(bad_site['fragment'])
            i,j,k,l = map(bad_frag.get_site_frag_idx, bad_site['site'])
            bad = bad_frag.solver.get_P2()[i,j,k,l]

            if not _is_zero(bad - good):
                kennyloggins.warning("Failed 2e- check: %s" % str(check))


    def _dmet_iter(s, pot):
        s._make_pots(pot)

        for frag in s.part.fragments:
            frag.solver.solve(frag.h+frag.hpot, frag.V+frag.Vpot, nelec=frag.ni)

        obj = s._calc_obj()
        return obj

    def _get_fragment_by_name(s,name):
        for frag in s.part.fragments:
            if frag.name == name:
                return frag
        else:
            kennyloggins.error("Unable to find fragment named %s. This is likely due to a bad partition file" % name) 
            import sys; sys.exit(1)

    def optimize(s, guess=None):
        if guess is not None: 
            guess = numpy.array(guess)
            obj = s._dmet_iter(guess)
            kennyloggins.info("Objective: %s" % str(obj))
            import sys
            sys.exit(0)
        
        guess = numpy.zeros(s.part.optlen())

        optpot = nr.nr(s._dmet_iter, guess)
        s._dmet_iter(optpot) # Just to make sure
        s._check()

        s.e = 0
        s.e1 = 0
        s.e2 = 0
        s.ec = 0
        for frag in s.part.fragments:
            kennyloggins.info("Final 1RDM for fragment %s:\n%s\n" % (frag.name, str(frag.solver.get_P())))
            ontop = numpy.array([frag.solver.get_P2()[i,i,i,i] for i in range(2*frag.ni)])
            kennyloggins.info("Final 2RDM (diagonal) for fragment %s:\n%s\n" % (frag.name, str(ontop)))
            kennyloggins.info("Final 2RDM (0101) frag %s:\n%s\n" % (frag.name, frag.solver.get_P2()[0,1,0,1]))
            s._calc_energy(frag,s.fock)
#            kennyloggins.info("Fragment %s core energy:    %f" % (frag.name, frag.core))
            kennyloggins.info("Fragment %s 1e- energy:    %f" % (frag.name, frag.e1))
            kennyloggins.info("Fragment %s 2e- total energy: %f" % (frag.name, frag.e2))
            kennyloggins.info("Fragment %s total energy:  %f" % (frag.name, frag.e))
            kennyloggins.info("Fragment %s 1e- center energy:    %f" % (frag.name, frag.e1c))
            kennyloggins.info("Fragment %s 2e- center energy: %f" % (frag.name, frag.e2c))
            kennyloggins.info("Fragment %s center energy: %f" % (frag.name, frag.ec))
            s.e += frag.e
            s.ec += frag.ec
            s.e1 += frag.e1
            s.e2 += frag.e2
        s.e /= len(s.part.fragments)
        s.ec /= len(s.part.fragments)
        s.e1 /= len(s.part.fragments)
        s.e2 /= len(s.part.fragments)

#        print "Dressed Hamiltonian???"
#        print s.part.fragments[0].solver.H_red
        #kennyloggins.debug("Final density matrix:\n"+str(s.solver.get_P()))



        return optpot

