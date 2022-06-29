# J.K. Ryan    
# 15 June 2018 

# import necessary stuff for code
import numpy as np
import sympy as sym


from numpy import *

import math
from math import *

#  Post-processing functions
import getkernelcoeff
import evalBDRYkernel

from evalBDRYkernel import evalBDRYkernel
from getkernelcoeff import getkernelcoeff


#######################################################################

def boundarypp(p,ell,RS,zEvalj,lambar):


    # Kernel coefficients -- ensures 2*RS+1 moments through polynomial reproduction
    # These are the B-spline weights
    cgam = np.zeros((2*RS+1))
    cgam = getkernelcoeff(ell,RS,lambar)

    
    # Get quadrature points and weights in order to evaluate the post-processing integrals:  
    # Approximation is a polynomial of degree p, kernel is a
    # polynomial of degree ell-1  ==> we need p+ell-1 = 2gpts-1, where n is the number of points.
    # Hence, gpts=(p+ell)/2.  If p+ell is odd, we want the first integer >= (p+ell)/2, hence the
    # ceiling function.
    gpts = int(ceil(0.5*(p+ell+1)))
    z = np.zeros((gpts))
    w = np.zeros((gpts))
    # Gauss-Legendre (default interval is [-1, 1])
    z, w = np.polynomial.legendre.leggauss(gpts)

    # Evaluate the Legendre polynomials at p+1 points per element                                  
    LPz = np.zeros((p+1,gpts))
    for i in arange(p+1):
        if i==0:
            LPz[i][:] = 1.0
        elif i ==1:
            LPz[i][:] = z
        else:
            LPz[i][:] = (2*i-1)/i*z[:]*LPz[i-1][:]-(i-1)/i*LPz[i-2][:]

    

    # Post-processor support is (xbar - kernelsupp*dx, xbar + kernelsupp*dx)
    kernelsupp = np.float(RS+0.5*ell)
    
    ksuppleft = int(ceil(-(kernelsupp+lambar)))
    
    # Make the element counter and integer value
    kwide = int(ceil(kernelsupp))
    # Total number of elements in the support
    pwide = 2*RS+ell
#    print('In boundarypp:\n')
#    print('    Evaluation point = ',zEvalj,'  Left support boundary=',ksuppleft,'\n')
#    print('    kernelsupp = ',kernelsupp,'    kwide = ',kwide,'    pwide = ',pwide,'\n')

    # bdrycc is the boundary post-processing matrix
    bdrycc = np.zeros((pwide,p+1))

    for kk1 in arange(pwide):
        kk = kk1+ksuppleft
        # Integral evaluation arrays
        ahat = np.float(-1.0)
        bhat = np.float(1.0)

        # Evaluation coordinate for the kernel integration                      
        kerzeta = np.zeros((gpts))
        
        kerzeta[:] = 0.5*(zEvalj-z[:])-np.float(kk)

        # Obtain the kernel value at the gauss points                           
        fker = np.zeros((gpts))
        fker = evalBDRYkernel(ell,RS,cgam,gpts,kerzeta,lambar)

        # Obtain the integral value     
        xintsum = np.zeros((p+1))                                        
        for m in arange(p+1):
            integralval = sum(fker[n]*LPz[m][n]*w[n] for n in arange(gpts))
            xintsum[m] = integralval


        # form the post-processing matrix = 0.5*(I1+I2)
        for m in arange(p+1):
            bdrycc[kk1][m] = 0.5*xintsum[m]

    return bdrycc


########################################################################  
