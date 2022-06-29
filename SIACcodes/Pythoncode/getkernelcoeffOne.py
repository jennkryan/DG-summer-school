# J.K. Ryan
# 18 May 2021

# import necessary stuff for code                                                                                                                                   
import numpy as np
import sympy as sym


from numpy import *
from scipy import *
from scipy import integrate
from scipy.special import binom
#  import matplotlib.pyplot as plt
import math
from math import *
#import np.linalg                                                                                                                                                   
import scipy.linalg   # SciPy Linear Algebra Library                                                                                                                
from scipy.linalg import lu
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve

import bspline

#######################################################################  

# Obtain the B-spline weights that give the kernel coefficients.  
# This is done through polynomial reproduction:
# int_R K(x-y)y^m dy = x^m, m=0..,2*RS.
# If the B-spline order is large, this matrix become ill-conditioned.

def getkernelcoeffOne(ell,RS,eta):
    # Define knot matrix
    
    print('In getkernelcoeffOne')
    width = 0
    
    tin = 1.e-14
    
    if abs(eta) > 0.5*(2*RS+ell)-tin:
    # symmetric kernel with 2*RS + 1 B-splines
        width = 0
    else:
        # boundary kernel with 2*RS + 2 B-splines
        width = 1
        
    RSOne = 2*RS+1+width
    
    # ng is the number of gauss points for quadrature = Bspline order * highest Number of moments
    ng = int(ceil(((ell-1)*RSOne-1)/2))
    xg = np.zeros((ng))
    wg = np.zeros((ng))
    
    print('eta = ',eta,'   ell = ',ell,'   RSone=',RSOne,'   ng=',ng,'\n')
    
    # Quadrature weights and nodes
    xg, wg = np.polynomial.legendre.leggauss(ng)
        
    
    T = np.zeros((RSOne,ell))
    for i in np.arange(RSOne):
        for j in np.arange(ell):
            if abs(eta)>0.5*(ell+2*RS)-tin:
            # nodes for symmetric filter
                T[i][j] = -(ell+RS)*0.5+j+i
            elif (eta<tin):
                # nodes for right boundary filter
                if i == 0 && j <= ell-1:
                    T[i][j] = eta
                elif  i == 1 && j== ell:
                    T[i][j] = eta + 1
                else:
                    T[i][j] = i+j-1+eta
            elif (eta > tin):
            # nodes for left boundary filter
                if i==RSOne:
                    if j==0
                        T[i][j] = eta-1
                    else:
                        T[i][j] = eta
                else
                    T[i][j] = -(ell+RS)+i+j+eta
            
    print('\n')
    print('Knot Matrix')
    print(T)
    print('\n')
     
    c = np.zeros((RSOne))
    A = np.zeros((RSOne,RSOne))
    for m in np.arange(RSOne):
        basis = bspline.Bspline(T[m][0:ell],ell-1)
        for j in np.arange(RSOne):
            component = 0
            for n in np.arange(m+1):
                jsum = 0.
                jsum = sum((-1)**(k+ell-1)*binom(ell-1,k)*((k-0.5*(ell-2))**(ell+n)-(k-0.5*ell)**(ell+n)) for k in np.arange(ell))
                component += binom(m,n)*(j-RS)**(m-n)*factorial(n)/factorial(n+ell)*jsum

            A[m][j] = component
            for k in np.arange(ell-1):
                xm = 0.5*(T[j][k]+T[j][k+1])
                xr = T[j][k+1] - T[j][k]
                for i in np.arange(ng):
                    A[m][k] = A[m][j] + xr*wg[i]*basis(xm+xr*xg[i])*(-(xr*xg[i]+xm))**(m-1)
    
    
    B = la.solve(A,B)
    
    b=np.zeros((RSOne))
    b[0]=1.

        
    #call the lu_factor function LU = linalg.lu_factor(A)
    Piv = scipy.linalg.lu_factor(A)
    #P, L, U = scipy.linalg.lu(A)
    #solve given LU and B
    c = scipy.linalg.lu_solve(Piv, b)
            

    print('\n')
    print('Matrix for SIAC coefficients:  \n')
    print(A)
    print('\n')

    print('SIAC coefficients:',c)

    # check coefficients add to one
    sumcoeff = sum(c[nel][n] for n in np.arange(2*RS+1))
    print('Sum of coefficients',sumcoeff)
    print('\n')
    print('\n')


    return c


########################################################################  
