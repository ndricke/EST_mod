import numpy
import scipy.optimize
import os
import sys
import logging
kennyloggins = logging.getLogger('nr')

import config

def fd_jacobian(f,x,eps=1e-6):
    n = len(x)
    J = numpy.zeros((n,n))
    for i in range(n):
        d = numpy.zeros(n)
        d[i] = 1
        ff  = f(x+d*eps)
        fff = f(x+2*d*eps)
        fb  = f(x-d*eps)
        fbb = f(x-2*d*eps)
        #J[:,i] = (f(x+d*eps) - f(x-d*eps))/(2*eps)
        #J[:,i] = (ff - fb)/(2*eps)
        J[:,i] = (-fff + 8*ff - 8*fb + fbb)/(12*eps)
    return J

def nr(f, x0):
    converged = False
    x = x0.copy()
    ite = 0

    while not converged:
        ite += 1
        fx = f(x)
        config.solvers['recycle_guess'] = True # Start using old points as guess
        J = fd_jacobian(f,x)

        #dx = numpy.linalg.solve(J,-fx)
        try:
            Jinv = numpy.linalg.inv(J)
        except(numpy.linalg.linalg.LinAlgError):
            print J
            kennyloggins.warning("Singular Jacobian:\n%s\n" % str(J))
            kennyloggins.warning("I'mma let you finish")
            return x
#            raise ValueError()
#            kennyloggins.warning("Attempting pseudoinverse")
#            Jinv = numpy.linalg.pinv(J)
        dx = Jinv.dot(-fx)
        
        #print "Jacobian"
        #print J
        #print 
        #print "Search direction:"
        #print dx
        #print

        scale = scipy.optimize.minimize_scalar(lambda s: numpy.linalg.norm(f(x+s*dx)),  tol=config.optimize['linesearch_tol'], method='bounded', bracket=[-1,2], bounds=[-1,2])
        if scale.x < 0:
            kennyloggins.warning("Stepped backwards")
        config.solvers['recycle_guess'] = False #Stop using old points as guess


        if numpy.linalg.norm(dx) < 1e-9 or numpy.linalg.norm(fx) < 1e-5:
            converged = True

        step = scale.x
        if config.optimize['max_step'] is not None:
            norm = numpy.linalg.norm(dx)
            if step*norm > config.optimize['max_step']:
                kennyloggins.warning("Limiting step size to max_step")
                step = config.optimize['max_step']/norm

        x += step * dx
        #x += dx
        kennyloggins.info("Iteration: " + str(ite))
        kennyloggins.debug("Jacobian\n" + 
            "---------------\n" + 
            str(J) + "\n" +
            str(numpy.linalg.inv(Jinv)))
        kennyloggins.info("x: " + str(x))
        kennyloggins.info("dx:" + str(dx) + " " + str(numpy.linalg.norm(dx)))
        kennyloggins.info("fx: " + str(fx) + " " + str(numpy.linalg.norm(fx)))
        kennyloggins.info("scale: " + str(scale.x))
        kennyloggins.info("scale iterations: " + str(scale.nfev))

    return x

    
