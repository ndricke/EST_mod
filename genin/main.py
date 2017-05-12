#! /usr/bin/env python

import argparse
import numpy
import logging
import sys
root = logging.getLogger()
root.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
root.addHandler(ch)

kennyloggins = logging.getLogger('main')

import hf
import lattice
import dmet
import solvers
import config
import partition

#numpy.set_printoptions(precision = 3, suppress = True)

parser = argparse.ArgumentParser()
parser.add_argument('-m', help='Number of electron pairs', type=int, required=True)
parser.add_argument('-p', '--partition',  help='Partition file', type=str, required=True)
parser.add_argument('-l', '--lattice',  help='Lattice file', type=str, required=True)
parser.add_argument('--core', help='Use core Hamiltonian correction', type=bool, default=True)
parser.add_argument('--solver', help="Impurity solver (GK|Troy)", type=str, default="Troy")
parser.add_argument('--debug', help="Debug output", action='store_true')
parser.add_argument('--guess', type=float, nargs='+')
parser.add_argument('--hfguess', type=str, default="None")
parser.add_argument('--no-two', help='Disable two electron checks and constraints', action='store_true')
args = parser.parse_args()

if args.debug:
    ch.setLevel(logging.DEBUG)

kennyloggins.debug("Input:\n%s" % args)

config.bootstrap['core'] = args.core
config.bootstrap['no_two'] = args.no_two

config.optimize['guess'] = args.guess
config.hf = args.hfguess

config.solvers['type'] = args.solver
l = lattice.Lattice.fromfile(args.lattice, args.m)
part = partition.Partition.fromfile(args.partition)

emb = None
kennyloggins.info("Bootstrap parameters:")
kennyloggins.info("Lattice: " + args.lattice)
kennyloggins.info("Filling: " + str(args.m))
kennyloggins.info("Partition: " + args.partition)
kennyloggins.info("Performing bootstrap optimization")
emb = dmet.Bootstrap(l, part)
emb.optimize(guess=args.guess)

kennyloggins.info("Final energy: %f" % emb.energy())
kennyloggins.info("Center energy: %f" % emb.center_energy())

kennyloggins.info("One electron energy: %f" % emb.e1)
kennyloggins.info("Two electron energy: %f" % emb.e2)

kennyloggins.info("Done")
logging.shutdown()

