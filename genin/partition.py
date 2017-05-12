import numpy
import os
import logging
import config
import collections

import fragment

kennyloggins = logging.getLogger('partition')

class Partition:
    @staticmethod
    def fromfile(fname):
        kennyloggins.debug("Loading partition from file " + os.path.abspath(fname))
        fragments = []

        constraints = {}
        constraints['filling'] = []
        constraints['1e'] = []
        constraints['2e'] = []
        constraints['2e_diff'] = []
        constraints['2e_energy'] = []
        constraints['2RDM_norm'] = []

        checks = {}
        checks['filling'] = []
        checks['1e'] = []
        checks['2e'] = []

        center = []

        f = open(fname)
        mode = None
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if line.startswith('#'):
                continue

            if not line:
                mode = None
            elif mode is None:
                if line.startswith("Fragment"):
                    mode = "Fragment"
                    fragments.append(fragment.Fragment(line.split()[1]))
                elif line == "Filling Constraints":
                    mode = "Filling Constraints"
                elif line == "1e- Constraints":
                    mode = "1e- Constraints"
                elif line == "2e- Constraints":
                    mode = "2e- Constraints"
                elif line == "2e- Difference Constraints":
                    mode = "2e- Difference Constraints"
                elif line == "2e- Energy Constraints":
                    mode = "2e- Energy Constraints"
                elif line == "Filling Checks":
                    mode = "Filling Checks"
                elif line == "1e- Checks":
                    mode = "1e- Checks"
                elif line == "2e- Checks":
                    mode = "2e- Checks"
                elif line == "2RDM Normalization Constraints":
                    mode = "2RDM Normalization Constraints"
                else:
                    kennyloggins.warning("Unrecongnized delimiter block: %s" % line)

            elif mode == "Fragment":
                last_frag = fragments[-1]
                if line.startswith("Center"):
                    last_frag.center = line.split()[1].split(',')
                else:
                    s = line.split()
                    last_frag.add_site(s[0], int(s[1]))

            elif mode == "Filling Constraints":
                sites_str = line.split(';')
                sites = map(Partition._parse_site, sites_str)
                constraints['filling'].append(sites)

            elif mode == "1e- Constraints":
                constraint = {}
                bads_str, good_str = line.split()
                constraint['good'] = Partition._parse_sites(good_str)
                constraint['bad'] = Partition._parse_many_sites(bads_str)
                constraints['1e'].append(constraint)

            elif mode == "2e- Constraints":
                constraint = {}
                bads_str, good_str = line.split()
                constraint['good'] = Partition._parse_sites(good_str)
                constraint['bad'] = Partition._parse_many_sites(bads_str)
                constraints['2e'].append(constraint)

            elif mode == "2e- Difference Constraints":
                constraint = {}
                plus_str, minus_str = line.split()
                constraint['plus'] = Partition._parse_many_sites(plus_str)
                constraint['minus'] = Partition._parse_many_sites(minus_str)
                constraints['2e_diff'].append(constraint)

            elif mode == "2e- Energy Constraints":
                constraint = {}
                bads_str, good_str = line.split()
                constraint['good'] = Partition._parse_site(good_str)
                constraint['bad'] = map(Partition._parse_site, bads_str.split(';'))
                constraints['2e_energy'].append(constraint)

            elif mode ==  "Filling Checks":
                checks['filling'].append(Partition._parse_site(line))

            elif mode == "1e- Checks":
                check = {}
                bad_str, good_str = line.split()
                check['good'] = Partition._parse_sites(good_str)
                check['bad'] = Partition._parse_sites(bad_str)
                checks['1e'].append(check)

            elif mode == "2e- Checks":
                check = {}
                bad_str, good_str = line.split()
                check['good'] = Partition._parse_sites(good_str)
                check['bad'] = Partition._parse_sites(bad_str)
                checks['2e'].append(check)

            elif mode == "2RDM Normalization Constraints":
                constraint = Partition._parse_many_sites(line)
                constraints['2RDM_norm'].append(constraint)
                print constraint

            else:
                kennyloggins.error("Unrecognized mode: %s\n This should never happen." % line)
                import sys; sys.exit(1)
        f.close()

        return Partition(fragments, constraints, checks)

    @staticmethod
    def _parse_site(string):
        s = string.split(':')
        site = {}
        site['fragment'] = s[0]
        site['site'] = s[1]
        return site

    @staticmethod
    def _parse_sites(string):
        s = string.split(':')
        site = {}
        site['fragment'] = s[0]
        site['site'] = s[1].split(',')
        return site

    @staticmethod
    def _parse_many_sites(string):
        return map(Partition._parse_sites, string.split(';'))


    def __init__(s, fragments, constraints, checks):
        s.fragments = fragments
        s.constraints = constraints
        s.checks = checks
        assert(len(s.fragments) >= 1)

        # check all fragments are uniquely named
        names = [frag.name for frag in s.fragments]
        if len(names) != len(set(names)):
            kennyloggins.error("Fragment name repeated. Error is likely in the partition file")
            import sys; sys.exit(1)

        # run indexing on each fragment
        for frag in s.fragments:
            frag.index()

        if '2RDM_norm' in s.constraints and len(s.constraints['2RDM_norm']) > 1:
            kennyloggins.error("Cannot have more that one 2RDM norm constraint")

        # remove disabled constraints and checks
        if config.bootstrap['no_two']:
            kennyloggins.info("Disabling 2e- constraints and checks")
            s.constraints['2e'] = []
            s.constraints['2e_diff'] = []
            s.constraints['2e_energy'] = []
            s.checks['2e'] = []
            s.constraints['2RDM_norm'] = []




    def optlen(s):
        return numpy.array([len(s.constraints[key]) for key in s.constraints]).sum()

