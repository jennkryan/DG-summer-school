# J.K. Ryan
# 18 May 2021

# import necessary stuff for code
import numpy as np
import sympy as sym

from numpy import *
from scipy import *
from scipy import integrate
from scipy.special import binom 

### work around for Mac OS Big Sur ###
import matplotlib as mpl
mpl.use('tkagg')
######################################
import matplotlib.pyplot as plt
import math
import scipy.linalg   # SciPy Linear Algebra Library

from scipy.linalg import lu
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve
from math import *

#  Post-processing functions
import getkernelcoeffOne

from getkernelcoeffOne import getkernelcoeffOne

########################################################################
# Test that the One-sided coefficients produced are correct
########################################################################  
#  MAIN PART OF PROGRAM
# Set number of evaluation points per element
evalPoints = 1

# Define number of elements, the polynomial degree
Nx = 5  #number of elements
#Nx = max(8,int(2*(2*p+1)))
#p = input('Input polynomial degree (NOT ORDER):  ')   #polynomial degree
#p = int(p)
p=1

# Get quadrature points and weights
z = np.zeros((p+1))
w = np.zeros((p+1))
# Gauss-Legendre (default interval is [-1, 1])
z, w = np.polynomial.legendre.leggauss(p+1)

# Now  get quadrature points and weights for the evaluation points
zEval = np.zeros((evalPoints))
wEval = np.zeros((evalPoints))
zEval, wEval = np.polynomial.legendre.leggauss(evalPoints)

# ASSUMING UNIFORM INTERVALS.  DOMAIN IS [0,1].
xright = np.float_(1.0)
xleft = np.float_(0.0)
xlength = xright - xleft
delta_x = xlength/float(Nx)
x_grid = np.zeros((Nx+1))
for k in np.arange(Nx+1):
    x_grid[k] = float(k)*delta_x

# Define kernel smoothness for post-processing
ellp2 = input('Input smoothness required (>=-1).  0 = continuous:  ');
ellp2 = int(ellp2)
# ell is the order of the B-spline
ell = ellp2 + 2

NumSpline = input('Input number of splines:  ');
NumSpline = int(NumSpline)
#RS = int(ceil((NumSpline-1)/2))
RS = ceil((NumSpline-1)/2)
print('RS = ',RS,'    Order of Spline = ',ell,'    p+ell=',p+ell,'\n')

# Half width of kernel support
kwide = int(ceil(RS+0.5*ell))
ksup = 2*kwide+1
print('kwide=',kwide,'    ksup=',ksup,'\n')


#symcc is the symmetric post-processing matrix
#symcc = symmetricpp(p,ell,RS,zEval)

# Calculate kernel coefficients based on location to boundary
# COne(element,eval pt)  is an array that stores the one-sided coefficients
# element = 0:(2*RS+ell)/2-1 is the left boundary
# element = (2*RS+ell)/2 : 2*RS+ell-1 is the right boundary

#COne = np.zeros((evalPoints,2*RS+ell+1,2*RS+2))
#etaLeft = np.zeros((evalPoints))
#etaRight = np.zeros((evalPoints))
for nel in np.arange(Nx):
    h = x_grid[nel+1] - x_grid[nel]
    zmap = 0.5*h*(zEval+1.0) + x_grid[nel]
    
    for j in np.arange(evalPoints):
        if nel < kwide:
        # Left boundary kernel coefficients
            etaLeft = (zmap[j]-xleft)/h
            #COne = getkernelcoeffOne(ell,RS,etaLeft)
            #left boundary kernel
            COne = getkernelcoeffOne(ell,RS,xleft,etaLeft)
            print('\n')
            print('nel = ',nel,'   Left Boundary Coefficients =')
            print(COne)
            print('\n')
        elif nel > Nx-2-kwide:
        # Right boundary kernel coefficients
            etaRight = (zmap[j]-xright)/h
            nelone = nel-(Nx-2-kwide)
            #COne[1:evalPoints][nelone] = getkernelcoeffOne(ell,RS,etaRight)
            # right boundary kernel
            COne = getkernelcoeffOne(ell,RS,xright,etaright)
            print('\n')
            print('nel = ',nel,'   Right Boundary Coefficients =')
            print(COne)
            print('\n')




