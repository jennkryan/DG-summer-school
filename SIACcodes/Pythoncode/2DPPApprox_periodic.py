# import necessary stuff for code
import numpy as np
import sympy as sym


from numpy import *
from scipy import *
from scipy import integrate
from scipy.special import binom 
import matplotlib.pyplot as plt
import math
#import np.linalg
import scipy.linalg   # SciPy Linear Algebra Library
from scipy.linalg import lu
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve

from math import *

import symmetricpp

from symmetricpp import symmetricpp


########################################################################


#  MAIN PART OF PROGRAM
# Set number of evaluation points per element
evalPoints = input('Input number of evaluation points per element:  ');
evalPoints = int(evalPoints)

# Define number of elements, the polynomial degree                                         
Nx = input('Input number of elements in x:  ');
Ny = input('Input number of elements in y:  ');
p = input('Input polynomial degree (NOT ORDER):  ');   #polynomial degree                  
Nx = int(Nx)
Ny = int(Ny)
p = int(p)

# Get quadrature points and weights

z = np.zeros((p+1))
w = np.zeros((p+1))
# Gauss-Legendre (default interval is [-1, 1])
z, w = np.polynomial.legendre.leggauss(p+1)

# Now  get quadrature points and weights for the evaluation points
zEval = np.zeros((evalPoints))
wEval = np.zeros((evalPoints))
zEval, wEval = np.polynomial.legendre.leggauss(evalPoints)

# Evaluate the Legendre polynomials at p+1 points per element
LegPolyAtz = np.zeros((p+1,p+1))
LegMass = np.zeros((p+1))
for i in arange(p+1):
    if i==0:
        LegPolyAtz[i][:] = 1.0
    elif i ==1:
        LegPolyAtz[i][:] = z
    else:
        LegPolyAtz[i][:] = (2*i-1)/i*z[:]*LegPolyAtz[i-1][:]-(i-1)/i*LegPolyAtz[i-2][:]
        
    LegMass[i] = np.dot(w,LegPolyAtz[i][:]*LegPolyAtz[i][:])

#Evaluate the Legendre polynomials at the evaluation points
LegPolyAtzEval = np.zeros((p+1,evalPoints))
for i in arange(p+1):
    if i==0:
        LegPolyAtzEval[i][:] = 1.0
    elif i ==1:
        LegPolyAtzEval[i][:] = zEval
    else:
        LegPolyAtzEval[i][:] = (2*i-1)/i*zEval[:]*LegPolyAtzEval[i-1][:]-(i-1)/i*LegPolyAtzEval[i-2][:]

# ASSUMING UNIFORM INTERVALS.  DOMAIN IS [0,1].                                            
xright = np.float_(1.0)
xleft = np.float_(0.0)
xlength = xright - xleft
hx = xlength/Nx
x_grid = np.zeros((Nx+1))
for k in arange(Nx+1):
    x_grid[k] = xleft + float(k)*hx

ytop = xright
ybottom = xleft
ylength = ytop - ybottom
hy = ylength/Ny
y_grid = x_grid


# Define modes of approximation by performing an L2-projection onto the space of piecewise Legendre polynomials                                                                   
uhat = np.zeros((Nx,Ny,p+1,p+1))
for nelx in arange(Nx):
    term1 = np.zeros((p+1))
    term3 = np.zeros((p+1))
    for mx in arange(p+1):
        for ix in arange(p+1):
            xrg = 0.5*hx*(z[ix]+1.0) + x_grid[nelx]
            term1[mx] = term1[mx] + np.sin(2.0*math.pi*xrg)*LegPolyAtz[mx][ix]/LegMass[mx]*w[ix]
            term3[mx] = term3[mx] + np.cos(2.0*math.pi*xrg)*LegPolyAtz[mx][ix]/LegMass[mx]*w[ix]

    for nely in arange(Ny):
        term2 = np.zeros((p+1))
        term4 = np.zeros((p+1))
        for my in arange(p+1):
            for jy in arange(p+1):
                yrg = 0.5*hy*(z[jy]+1.0) + y_grid[nely]
                term2[my] = term2[my] + np.cos(2.0*math.pi*yrg)*LegPolyAtz[my][jy]/LegMass[my]*w[jy]
                term4[my] = term4[my] + np.sin(2.0*math.pi*yrg)*LegPolyAtz[my][jy]/LegMass[my]*w[jy]


        for mx in arange(p+1):
            for my in arange(p+1):
                uhat[nelx][nely][mx][my] = term1[mx]*term2[my] + term3[mx]*term4[my]




# Gather values in order to plot the exact solution and the projection                                                                                                            
xEval=[]
yEval=[]
fExact=[]
fApprox=[]

NxEval = Nx*evalPoints
NyEval = Ny*evalPoints


for nelx in arange(Nx):
    for ix in arange(evalPoints):

        xrg = 0.5*hx*(zEval[ix]+1.0) + x_grid[nelx]

        for nely in arange(Ny):
            for jy in arange(evalPoints):


                yrg = 0.5*hy*(zEval[jy]+1.0) + y_grid[nely]
                xEval.append(xrg)
                yEval.append(yrg)

                fval = np.sin(2.0*math.pi*(xrg+yrg))
                fExact.append(fval)

                uval = np.float_(0.0)
                for mx in arange(p+1):
                    for my in arange(p+1):
                        uval = uval + uhat[nelx][nely][mx][my]*LegPolyAtzEval[mx][ix]*LegPolyAtzEval[my][jy]


                fApprox.append(uval)



ApproxErr = np.subtract(fExact, fApprox)
PtwiseAErr = np.abs(ApproxErr)
LinfAErr = np.max(PtwiseAErr)


#############################################################################################################

# Define kernel smoothness
ellp2 = input('Input smoothness required (>=0).  0 = continuous:  ');
ellp2 = int(ellp2)
ell = ellp2 + 2

problemtype = input('Inpute 1=elliptic/parabolic, 2=hyperbolic:  ')
problemtype = int(problemtype)
# print(problemtype)
# If Elliptic/Parabolic -- Define the number of splines (2*RS+2)
if problemtype == 1:
    if p+ell >= 2*p+2:
        RS = p+1
    elif p+ell <=p+1:
        RS = int(ceil(p/2))+1
    else:
        RS = int(ceil((p+ell-1)/2))+1
elif problemtype == 2:
    # If Hyperbolic -- Define the number of splines (2*RS+1)
    # Define number of splines
    if p+ell >= 2*p+1:
        RS = p
    elif p+ell <=p+1:
        RS = int(ceil(p/2))
    else:
        RS = int(ceil((p+ell-1)/2))

kwide = ceil(RS+0.5*ell)

#symcc is the symmetric post-processing matrix                                 
symcc = symmetricpp(p,ell,RS,zEval)


PPxEval=[]
PPyEval=[]
PPfExact=[]
PPfApprox=[]
onedcoord=[]
f1d=[]
pp1d=[]


for nelx in arange(Nx):
    for ix in arange(evalPoints):
        xrg = 0.5*hx*(zEval[ix]+1.0) + x_grid[nelx]
        PPxEval.append(xrg)
        for nely in arange(Ny):
            for jy in arange(evalPoints):
                yrg = 0.5*hy*(zEval[jy]+1.0) + y_grid[nely]
                
                f=np.sin(2.0*math.pi*(xrg+yrg))
                PPyEval.append(yrg)
                PPfExact.append(f)

                # set indices to form post-processed solution
                upost = 0.0
                if kwide <= nelx <= Nx-2-kwide:
                    for kkx in arange(2*kwide+1):
                        kk2x = kkx - kwide
                        xindex = nelx + kk2x

                        if kwide <= nely <= Ny-2-kwide:
                            for kky in arange(2*kwide+1):
                                kk2y = kky - kwide
                                yindex = nely + kk2y
                                for mx in arange(p+1):
                                    for my in arange(p+1):
                                        upost = upost + symcc[kkx][mx][ix]*symcc[kky][my][jy]*uhat[xindex][yindex][mx][my]
        
                        elif nely < kwide:
                            for kky in arange(2*kwide+1):
                                kk2y = kky - kwide
                                if nely+kk2y<0:
                                    yindex = Ny+nely+kk2y
                                else:
                                    yindex = nely+kk2y
                                for mx in arange(p+1):
                                    for my in arange(p+1):
                                        upost = upost + symcc[kkx][mx][ix]*symcc[kky][my][jy]*uhat[xindex][yindex][mx][my]
                                        
                        else:
                            for kky in arange(2*kwide+1):
                                kk2y = kky - kwide
                                if kk2y <=0:
                                    yindex = nely + kk2y
                                else:
                                    yindex = nely-Ny+kk2y
                                for mx in arange(p+1):
                                    for my in arange(p+1):
                                        upost = upost + symcc[kkx][mx][ix]*symcc[kky][my][jy]*uhat[xindex][yindex][mx][my]
                                
                elif nelx < kwide:
                    for kkx in arange(2*kwide+1):
                        kk2x = kkx - kwide
                        if nelx+kk2x <0:
                            xindex = Nx+nelx+kk2x
                        else:
                            xindex = nelx+kk2x
                        if kwide <=nely <= Ny-2-kwide:
                            for kky in arange(2*kwide+1):
                                kk2y = kky - kwide
                                yindex = nely+kk2y
                                for mx in arange(p+1):
                                    for my in arange(p+1):
                                        upost = upost + symcc[kkx][mx][ix]*symcc[kky][my][jy]*uhat[xindex][yindex][mx][my]
                        elif nely < kwide:
                            for kky in arange(2*kwide+1):
                                kk2y = kky - kwide
                                if nely + kk2y < 0:
                                    yindex = Ny + nely + kk2y
                                else:
                                    yindex = nely + kk2y
                                for mx in arange(p+1):
                                    for my in arange(p+1):
                                        upost = upost + symcc[kkx][mx][ix]*symcc[kky][my][jy]*uhat[xindex][yindex][mx][my]
                        else:
                            for kky in arange(2*kwide+1):
                                kk2y = kky - kwide
                                if kk2y <=0:
                                    yindex = nely + kk2y
                                else:
                                    yindex = nely - Ny +kk2y
                                for mx in arange(p+1):
                                    for my in arange(p+1):
                                        upost = upost + symcc[kkx][mx][ix]*symcc[kky][my][jy]*uhat[xindex][yindex][mx][my]
                else:
                    for kkx in arange(2*kwide+1):
                        kk2x = kkx - kwide
                        if kk2x <= 0:
                            xindex = nelx + kk2x
                        else:
                            xindex = nelx-Nx+kk2x
                        if kwide <=nely <= Ny-2-kwide:
                            for kky in arange(2*kwide+1):
                                kk2y = kky - kwide
                                yindex = nely+kk2y
                                for mx in arange(p+1):
                                    for my in arange(p+1):
                                        upost = upost + symcc[kkx][mx][ix]*symcc[kky][my][jy]*uhat[xindex][yindex][mx][my]
                        elif nely < kwide:
                            for kky in arange(2*kwide+1):
                                kk2y = kky - kwide
                                if nely + kk2y < 0:
                                    yindex = Ny + nely + kk2y
                                else:
                                    yindex = nely + kk2y
                                for mx in arange(p+1):
                                    for my in arange(p+1):
                                        upost = upost + symcc[kkx][mx][ix]*symcc[kky][my][jy]*uhat[xindex][yindex][mx][my]
                        else:
                            for kky in arange(2*kwide+1):
                                kk2y = kky - kwide
                                if kk2y<=0:
                                    yindex = nely + kk2y
                                else:
                                    yindex = nely - Ny +kk2y
                                for mx in arange(p+1):
                                    for my in arange(p+1):
                                        upost = upost + symcc[kkx][mx][ix]*symcc[kky][my][jy]*uhat[xindex][yindex][mx][my]

                                     
                PPfApprox.append(upost)
#                if (nelx == 4) and (ix == 0):
#                    onedcoord.append(yrg)
#                    f1d.append(f)
#                    pp1d.append(upost)


#ppdiff = np.subtract(f1d,pp1d)
#ppSliceerr = np.abs(ppdiff)


PPApproxErr = np.subtract(PPfExact, PPfApprox)
PPPtwiseAErr = np.abs(PPApproxErr)
PPLinfAErr = np.max(PPPtwiseAErr)


print('\n')
print('Nx = ',Nx,'    Ny = ',Ny,'    p = ',p)
print('L-inf Error for the Projection =',LinfAErr)
print('L-inf Error for the Post-processed Projection =',PPLinfAErr)
print('\n')


#plt.figure(1)
#plt.plot(onedcoord,f1d,y1D,Approxslice,onedcoord,pp1d)
#plt.legend(['Exact','L2-Projection','SIAC DG'])

#plt.figure(2)
#plt.semilogy(y1D,Sliceerr,onedcoord,ppSliceerr)

#plt.show()



