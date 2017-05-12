import logging
import config
import solvers
kennyloggins = logging.getLogger('fragment')

class Fragment:
    def __init__(s, name):
        s.name = name
        if config.solvers['type'].lower() == "gk":
            s.solver = solvers.GKFCI()
        elif config.solvers['type'].lower() == "troy":
            s.solver = solvers.TroyFCI()
        else:
            kennyloggins.error("Unknown solver: %s" % config.solvers['type'])
            import sys; sys.exit(1)
        s.sites = []
        s.center = []
        s.indexed = False

    def add_site(s, name, lat_idx):
        site = {}
        site['index'] = int(lat_idx)
        site['name'] = name
        s.sites.append(site)


    def index(s):
        s.name_to_frag_idx = {item['name']:i for i,item in enumerate(s.sites)}
        s.name_to_lat_idx = {item['name']:item['index'] for i,item in enumerate(s.sites)}
        if len(s.name_to_frag_idx) != len(s.sites):
            kennyloggins.error("Site name repeated. Error is likely in the partition file.")
            import sys; sys.exit(1)
        s.lat_idx = [site['index'] for site in s.sites]
        s.c_lat_idx = [s.name_to_lat_idx[ci] for ci in s.center]
        s.ni = len(s.sites)

        s.indexed = True

    def get_site_frag_idx(s, name):
        if not s.indexed:
            kennyloggins.error("Fragment %s not indexed!" % s.name)
            import sys; sys.exit(1)
        return s.name_to_frag_idx[name]

    def get_site_lat_idx(s, name):
        if not s.indexed:
            kennyloggins.error("Fragment %s not indexed!" % s.name)
            import sys; sys.exit(1)
        return s.name_to_lat_idx[name]

    def lattice_indices(s):
        if not s.indexed:
            kennyloggins.error("Fragment %s not indexed!" % s.name)
            import sys; sys.exit(1)
        return s.lat_idx

    def center_lattice_indices(s):
        if not s.indexed:
            kennyloggins.error("Fragment %s not indexed!" % s.name)
            import sys; sys.exit(1)
        return s.c_lat_idx

