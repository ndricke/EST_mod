import numpy
import scipy.optimize
import os
import sys
import logging
kennyloggins = logging.getLogger('al')

ctol = 1e-5
mu_multiplier = 5
mu_init = 10

def objective(f, x, mu, la):
    fx = f(x)
    return fx[0] + mu/2*(fx[1]**2).sum() - (la*fx[1]).sum()

def al(f, x0):
    kennyloggins.info("Beginning augmented lagrangian optimization")
    kennyloggins.debug("Constraint tolerence: %f" % ctol)
    converged = False
    x = x0.copy()
    ite = 0
    
    mu = mu_init
    fx = f(x) # XXX: determining number of constraints by length of output
    la = numpy.zeros(len(fx[1]))
    kennyloggins.debug("Initial mu: %f" % mu)
    kennyloggins.debug("Initial la: " + str(la))


    # outer loop: update mu and la
    while not converged:
        ite += 1
        kennyloggins.info("Outer Iteration %d" % ite)
        kennyloggins.debug("Mu = %f" %mu)
        kennyloggins.debug("La = %s" % str(la))

        x = scipy.optimize.fmin_bfgs(lambda xi: objective(f,xi,mu=mu,la=la), x)
        kennyloggins.debug("Optimized x: " + str(x))
        fx = f(x)
        kennyloggins.info("Objective: %f" % fx[0])
        kennyloggins.info("Worst constraint: %f" % numpy.abs(fx[1]).max())

        if numpy.abs(fx[1]).max() < ctol:
            converged = True
            kennyloggins.info("Converged")
            break
        la -= mu*fx[1]
        mu *= mu_multiplier
    kennyloggins.debug("Final mu: %f" % mu)
    kennyloggins.debug("Final la: " + str(la))

    kennyloggins.info("Final x: " + str(x))

    return x

    
