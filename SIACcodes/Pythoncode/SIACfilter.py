# import necessary stuff for code
# change
import numpy as np
import sympy as sym

from numpy import *
from scipy import *
from scipy import integrate
from scipy.special import binom

import os
import math
#import np.linalg
import scipy.linalg   # SciPy Linear Algebra Library
from scipy.linalg import lu
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve

from math import *

import bspline
from DG_Approximation import DGScheme

########################################################################

#class SIACfilter(siac_config,smoothness,polynomial_degree,num_grid_cells):
class SIACfilter(object):
    # Pass in order of dg approximation and smoothness of reconstruction.
    # order>0 and smoothness>=0
    # The reconstruction is computed at the points given by zEval
    # in the reference element [-1,1]
    def __init__(self,grid,**kwargs):
        # Unpack keyword arguments
        self._grid = grid
        self._polynomial_degree = kwargs.pop('polynomial_degree', 2)
        self._num_grid_cells = kwargs.pop('num_grid_cells', 16)
        self._left_bound = kwargs.pop('left_bound', -1)
        self._right_bound = kwargs.pop('right_bound', 1)
        self._smoothness = kwargs.pop('smoothness',1)
        self._num_eval_points = kwargs.pop('num_eval_points', 1)
        self._ell = self._smoothness + 2
        #self._RS = int(max(ceil(0.5*(self.order_+self.ell_-1)),ceil(0.5*self.order_)));
        self._RS = max(1,self._polynomial_degree)
        self._kwide = int(ceil(self._RS+0.5*self._ell))
        #self._num_eval_points = len(zEval)
        #self._symcc = symmetricpp(self._polynomial_degree+1,self._ell,smoothness,self._RS,self._num_eval_points,self._zEval)
        #self._symcc = symmetricpp(self)

        self._num_quad_points = self._polynomial_degree + 1
        
        #self._L = L
        
        # Throw an error if there are extra keyword arguments
        if len(kwargs) > 0:
            extra = ', '.join('"%s"' % k for k in list(kwargs.keys()))
            raise ValueError('Unrecognized arguments: %s' % extra)
            
        self._reset()
            


    def get_name(self):
        return self.__class__.__name__

    # Evaluate the reconstruction at all points given in the constructor
    # for the grid cell 'nel'. uhat has to contain data for the cells
    # 'nel-kwide' to 'nel+kwide' otherwise an exception is raised.
    def evaluate(self,dg_scheme,grid,uhat):
    
        PPxEval = grid
        PPfApprox = np.zeros((self._num_eval_points*self._num_grid_cells,1))
    
        zEval,symcc = self.symmetricpp()
        
        print('symcc = ',symcc)
        for nel in arange(self._num_grid_cells):
            
            #if ( nel - self._kwide < 0 or nel + self._kwide >= len(uhat) ):
            #    raise ValueError("element stencil not large enough")
   
            upost = np.zeros((self._num_eval_points))
            if self._kwide <= nel <= self._num_grid_cells-2-self._kwide:
                # Post-process interior elements
                upost = sum(sum(symcc[kk][m][:]*uhat[m][nel+kk-self._kwide] for m in arange(self._polynomial_degree+1)) for kk in arange(2*self._kwide+1))
            elif nel < self._kwide:
                print('nel = ',nel,' left boundary')
                # Left boundary elements
                for j in arange(self._num_eval_points):
                    for kk in arange(2*self._kwide+1):
                        kk2 = kk-self._kwide
                        for m in arange(self._polynomial_degree+1):
                            # periodic
                            if nel+kk2 < 0:
                                if self._num_eval_points == 1:
                                    upost = upost + symcc[kk][m][0]*uhat[m][self._num_grid_cells+nel+kk2]
                                else:
                                    upost[j] = upost[j] + symcc[kk][m][j]*uhat[m][self._num_grid_cells+nel+kk2]
                            else:
                                if self._num_eval_points == 1:
                                    upost = upost + symcc[kk][m][0]*uhat[m][nel+kk2]
                                else:
                                    upost[j] = upost[j] + symcc[kk][m][j]*uhat[m][nel+kk2]
            elif nel > self._num_grid_cells-2-self._kwide:
                # Right boundary elements
                print('nel = ',nel,' right boundary')
                for j in arange(self._num_eval_points+1):
                    for kk in arange(2*self._kwide+1):
                        kk2 = kk-self._kwide
                        for m in arange(self._polynomial_degree+1):
                        # periodic
                            if kk2 <=0:
                                if self._num_eval_points == 1:
                                    upost = upost + symcc[kk][m]*uhat[m][nel+kk2]
                                else:
                                    upost[j] = upost[j]+symcc[kk][m][j]*uhat[m][nel+kk2]
                            else:
                                if self._num_eval_points == 1:
                                    upost = upost + symcc[kk][m][0]*uhat[m][nel-self._num_grid_cells+kk2]
                                else:
                                    upost[j] = upost[j]+symcc[kk][m][j]*uhat[m][nel-self._num_grid_cells+kk2]

            
            if self._num_eval_points == 1:
                PPfApprox[nel] = upost
            else:
                upost = np.transpose(np.array(upost))
                PPfApprox[nel*self._num_eval_points:(nel+1)*self._num_eval_points] = upost
                        
        return PPfApprox
        
    def symmetricpp(self):

        # Kernel coefficients -- ensures 2*RS+1 moments through polynomial reproduction
        # These are the B-spline weights
        cgam = self.getkernelcoeff()

    
        # Get quadrature points and weights in order to evaluate the post-processing integrals:
        # Approximation is a polynomial of degree p, kernel is a
        # polynomial of degree ell-1  ==> we need p+ell-1 = 2gpts-1, where n is the number of points.
        # Hence, gpts=(p+ell)/2.  If p+ell is odd, we want the first integer >= (p+ell)/2, hence the
        # ceiling function.
        gpts = int(ceil(0.5*(self._polynomial_degree+self._ell)))
        z = np.zeros((gpts))
        w = np.zeros((gpts))
        # Gauss-Legendre (default interval is [-1, 1])
        z, w = np.polynomial.legendre.leggauss(gpts)
        
        zEval = np.zeros((self._num_eval_points))
        wEval = np.zeros((self._num_eval_points))
        zEval, wEval = np.polynomial.legendre.leggauss(self._num_eval_points)

        # Post-processor support is (xbar - kernelsupp*dx, xbar + kernelsupp*dx)
        kernelsupp = np.float(self._RS+0.5*self._ell)
        # Make the element counter and integer value
        #kwide = int(ceil(kernelsupp))
        # Total number of elements in the support
        pwide = 2*self._kwide+1
        print('In Symmetricpp:\n')
        print('    kernelsupp = ',kernelsupp,'    kwide = ',self._kwide,'    pwide = ',pwide,'\n')

        # Need to account for the case where the B-spline breaks are not aligned with the evaluation point.
        # This occurs with ell is odd (odd B-spline order)
        if self._ell % 2 == 0:
            kres = np.float(0)
        else:
            kres = np.float(1.0)

        #   symcc is the symmetric post-processing matrix
        symcc = np.zeros((pwide,self._polynomial_degree+1,self._num_eval_points))

        for j in arange(self._num_eval_points): # Loop over element evaluation points

            if kres !=0 and zEval[j] > 0: # locate kernel break.  Done based on where the
                                      # evaluation point is  with respect to cell center
                                      # if ell is odd
                kres = np.float(-1.0)
        
            zetaEval = zEval[j]+kres # This is the location of the kernel break within the element
                                 # for a uniform grid
        
            for kk1 in arange(pwide):
                kk = kk1 - self._kwide
            # Integral evaluation arrays
        
                if self._ell % 2 == 0: #B-spline order is even.
                             #Kernel breaks are the shifted evaluation point

                    # evaluation od the first integral
                    # in_{-1}^zEval[j] K(0.5*(zEval[j]-x)-kk)P[m][x] dx, where kk is
                    # the current element with respect the element of the post-processing point
                    ahat = np.float(-1.0)
                    bhat = np.float(zetaEval)

                    xintsum1 = self.xintsum(ahat,bhat,pwide,kernelsupp,z,w,cgam,kk,zetaEval,zEval[j])
                    # evaluation of the second integral (later need to scale by 2)
                    # int_zEval[j]^1 K(0.5*(zEval[j]-x)-kk)P[m][x] dx, where
                    # kk is the current element with respect the element of
                    # the post-processing point
                    ahat = np.float(zetaEval)
                    bhat = np.float(1.0)

                    xintsum2 =  self.xintsum(ahat,bhat,pwide,kernelsupp,z,w,cgam,kk,zetaEval,zEval[j])
                else:  #B-spline order is odd, kernel breaks depend on location of evaluation point
                   #with respect to the cell center
                
                    if zEval[j] != 0:
                        ahat = np.float(-1.0)
                        bhat = np.float(zetaEval)
                        xintsum1 =  self.xintsum(ahat,bhat,pwide,kernelsupp,z,w,cgam,kk,zetaEval,zEval[j])
                        #xintsum1 = xintsum(ahat,bhat,p,kwide,pwide,kernelsupp,z,w,RS,ell,cgam,kk,zetaEval,zEval[j])

                        ahat = np.float(zetaEval)
                        bhat = np.float(1.0)
                        xintsum2 =  self.xintsum(ahat,bhat,pwide,kernelsupp,z,w,cgam,kk,zetaEval,zEval[j])

                    else:
                        ahat = np.float(-1.0)
                        bhat = np.float(1.0)
                        xintsum1 =  self.xintsum(ahat,bhat,pwide,kernelsupp,z,w,cgam,kk,zetaEval,zEval[j])
                        xintsum2 = np.zeros((p+1))

                # form the post-processing matrix = 0.5*(I1+I2)
                for m in arange(self._polynomial_degree+1):
                    symcc[kk1][m][j] = 0.5*(xintsum1[m]+xintsum2[m])

        return zEval, symcc

    # Obtain the B-spline weights that give the kernel coefficients.
    # This is done through polynomial reproduction:
    # int_R K(x-y)y^m dy = x^m, m=0..,2*RS.
    # If the B-spline order is large, this matrix become ill-conditioned.

    def getkernelcoeff(self):
        # Define matrix to determine kernel coefficients
        A=np.zeros((2*self._RS+1,2*self._RS+1))
        for m in arange(2*self._RS+1):
            for gam in arange(2*self._RS+1):
                component = 0.
                for n in arange(m+1):
                    jsum = 0.
                    jsum = sum((-1)**(j+self._ell-1)*binom(self._ell-1,j)*((j-0.5*(self._ell-2))**(self._ell+n)-(j-0.5*self._ell)**(self._ell+n)) for j in arange(self._ell))
                    component += binom(m,n)*(gam-self._RS)**(m-n)*factorial(n)/factorial(n+self._ell)*jsum

                    A[m][gam] = component

        print('\n')
        print('Matrix for SIAC coefficients')
        print(A)
        print('\n')

        b=np.zeros((2*self._RS+1))
        b[0]=1.

        c = np.zeros((2*self._RS+1))
        #call the lu_factor function LU = linalg.lu_factor(A)
        Piv = scipy.linalg.lu_factor(A)
        #P, L, U = scipy.linalg.lu(A)
        #solve given LU and B
        c = scipy.linalg.lu_solve(Piv, b)


        print('SIAC coefficients:',c)

        # check coefficients add to one
        sumcoeff = sum(c[n] for n in arange(2*self._RS+1))
        print('Sum of coefficients',sumcoeff)


        return c
        
        
    # Evaluate the post-processing integrals using Gauss-Legendre quadrature.  The integral is:
    # int_a^b K(0.5(zEval - x) - kk)P^(m)(x) dx,
    # where K is the SIAC kernel using 2*RS+1 Bsplines of order ell,
    # zEval is the evaluation point, and P^(m) is the Legendre polynomial of degree m.

    def xintsum(self,ahat,bhat,pwide,kernelsupp,z,w,cgam,kk,zetaEval,point):

        gpts = int(len(z))
        xintsum = np.zeros((self._polynomial_degree+1))

        # Ensure integration does not go beyond the support of the kernel.
        intlow = np.float(point-2*kernelsupp-2.0*kk)
        intup = np.float(point+2*kernelsupp-2.0*kk)
        if ahat < intlow:
            ahat = intlow
        if bhat > intup:

            bhat = intup

        if ahat < bhat: # only perform the integration if ahat < bhat
                    # scale the integration interval to (-1,1) to use Gauss-Legendre quadrature
            abplus = 0.5*(ahat+bhat)
            abminus = 0.5*(bhat-ahat)
            zeta = np.zeros((gpts))
            zeta[:] = abminus*z[:]+abplus # quadrature coordinate

            # Evaluation coordinate for the kernel integration
            kerzeta = np.zeros((gpts))
            kerzeta[:] = 0.5*(point-zeta[:])-np.float(kk)

            # Obtain the kernel value at the gauss points
            fker = np.zeros((gpts))
            fker = self.evalkernel(cgam,gpts,kerzeta)

            # Legendre polynomials evaluated at the gauss points
            PLeg = np.zeros((self._polynomial_degree+1,gpts))
            for m in arange(self._polynomial_degree+1):
                if m==0:
                    PLeg[m][:] = np.ones((gpts))
                elif m ==1:
                    PLeg[m][:] = zeta
                else:
                    PLeg[m][:] = (2.0*m-1.0)/np.float(m)*zeta[:]*PLeg[m-1][:]-(m-1.0)/np.float(m)*PLeg[m-2][:]
            
            # Obtain the integral value
            for m in arange(self._polynomial_degree+1):
                integralval = sum(fker[n]*PLeg[m][n]*w[n] for n in arange(gpts))
                xintsum[m] = abminus*integralval
        
        return xintsum
    
                               
    #  evaluate the kernel:
    #  K(x) = sum_gam c_gam psi^{ell}(x-gam),
    #  where psi^{ell}(x-gam) is a B-spline of order ell centered at gam.

    def evalkernel(self,cgam,gpts,kerzeta):

        # Define B-spline breaks for a B-spline of order ell
        bsbrks = np.linspace(-0.5*(2*self._RS+self._ell),0.5*(2*self._RS+self._ell),2*self._RS+self._ell+1)
        basis = bspline.Bspline(bsbrks,self._ell-1)

    #    basis.plot()
    #    bsbrks = np.zeros(ell+1)
    #    for i in arange(ell+1):
    #        bsbrks[i] = -0.5*ell+i


        fker = np.zeros((gpts))

        g = np.zeros((gpts, 2*self._RS+1))

        for n in arange(gpts): # summing over zetas
            g[n][:] = basis(kerzeta[n])*cgam[:]
    
            fker[n] = sum(g[n][jj] for jj in arange(2*self._RS+1))
                
        return fker

########################################################################
    def _reset(self):
        # Set additional necessary instance variables
        #self._right_bound = self._right_bound.astype('float')
        #self._left_bound = self._left_bound.astype('float')
        self._interval_len = self._right_bound-self._left_bound
        #self._interval_len = self._interval_len.astype('float')
        self._cell_len = self._interval_len / self._num_grid_cells
        #self._basis = OrthonormalLegendre(self._polynomial_degree)

        # Set mesh with one ghost point on each side
        self._mesh = np.arange(self._left_bound - (3/2*self._cell_len), self._right_bound + (5/2*self._cell_len), self._cell_len)  # +3/2
