# import necessary stuff for code
# change
import numpy as np
import sympy as sym
import os
import math

#from numpy import *
#from scipy import *
from scipy import integrate
from scipy.special import binom


import matplotlib.pyplot as plt
import seaborn as sns
#import math
#import np.linalg
#import scipy.linalg   # SciPy Linear Algebra Library
#from scipy.linalg import lu
#from scipy.linalg import lu_factor
#from scipy.linalg import lu_solve

#from math import *

import Initial_Condition
from Basis_Function import OrthonormalLegendre
#import bspline
#from DG_Approximation import DGScheme
#from SIACfilter import SIACfilter
sns.set_palette("colorblind")


x = sym.Symbol('x')
z = sym.Symbol('z')

########################################################################

#class SIACfilter(siac_config,smoothness,polynomial_degree,num_grid_cells):
class SIACplot(object):
#(init_cond,init_config,polynomial_degree,num_grid_cells,left_bound,right_bound,num_eval_points,smoothness):
    # Pass in order of dg approximation and smoothness of reconstruction.
    # order>0 and smoothness>=0
    # The reconstruction is computed at the points given by zEval
    # in the reference element [-1,1]
    def __init__(self,**kwargs):
        # Unpack keyword arguments
        self._init_cond = kwargs.pop('init_cond', 'Sine')
        self._init_config = kwargs.pop('init_config', {})
        self._polynomial_degree = kwargs.pop('polynomial_degree', 2)
        self._num_grid_cells = kwargs.pop('num_grid_cells', 16)
        self._left_bound = kwargs.pop('left_bound', -1)
        self._right_bound = kwargs.pop('right_bound', 1)
        #self._smoothness = kwargs.pop('smoothness',1)
        self._num_eval_points = kwargs.pop('num_eval_points', 1)
        self._smoothness = kwargs.pop('smoothness',1)
        self._RS = max(1,self._polynomial_degree)
        
        
        #self._ell = self._smoothness + 2
        #self._RS = int(max(ceil(0.5*(self.order_+self.ell_-1)),ceil(0.5*self.order_)));
        #self._RS = max(1,self._polynomial_degree)
        #self._kwide = int(math.ceil(self._RS+0.5*self._ell))
        #self._num_eval_points = len(zEval)
        #self._symcc = symmetricpp(self._polynomial_degree+1,self._ell,smoothness,self._RS,self._num_eval_points,self._zEval)
        #self._symcc = symmetricpp(self)
        
        #self._num_quad_points = self._polynomial_degree + 1
        
        #self._L = L
        # Set parameters from config if existing
        self._plot_dir = kwargs.pop('plot_dir', 'fig')
        self._colors = kwargs.pop('colors', {})
        
        # Throw an error if there are extra keyword arguments
        if len(kwargs) > 0:
            extra = ', '.join('"%s"' % k for k in list(kwargs.keys()))
            raise ValueError('Unrecognized arguments: %s' % extra)
            
        # Make sure all classes actually exist
        if not hasattr(Initial_Condition, self._init_cond):
            raise ValueError('Invalid initial condition: "%s"' % self._init_cond)
            
        self._check_colors()
        self._reset()
                
        # Replace the string names with the actual class instances
        self._init_cond = getattr(Initial_Condition, self._init_cond)(self._left_bound, self._right_bound, self._init_config)
        
    def _check_colors(self):
        self._colors['exact'] = self._colors.get('exact', 'k-')
        self._colors['approx'] = self._colors.get('approx','b:')
        self._colors['filtered'] = self._colors.get('filtered', 'r--')
        
    def _reset(self, config):
        sns.set()

    def get_name(self):
        return self.__class__.__name__
        
    def dg_save_plots(self):
        name = self._init_cond.get_name() + '__' +  '__number_of_cells_' \
            + str(self._num_grid_cells) + '__polynomial_degree_' + str(self._polynomial_degree)

        # Set paths for plot files if not existing already
        if not os.path.exists(self._plot_dir):
            os.makedirs(self._plot_dir)

        if not os.path.exists(self._plot_dir + '/exact_and_approx'):
            os.makedirs(self._plot_dir + '/exact_and_approx')

        if not os.path.exists(self._plot_dir + '/semilog_error'):
            os.makedirs(self._plot_dir + '/semilog_error')

        if not os.path.exists(self._plot_dir + '/error'):
            os.makedirs(self._plot_dir + '/error')

        # Save plots
        plt.figure(1)
        plt.savefig(self._plot_dir + '/exact_and_approx/' + name + '.pdf')

        plt.figure(2)
        plt.savefig(self._plot_dir + '/semilog_error/' + name + '.pdf')

        plt.figure(3)
        plt.savefig(self._plot_dir + '/error/' + name + '.pdf')

    def dg_plot_results(self, zEval,wEval,grid, projection):
        
        max_error = self._dg_plot_mesh(zEval,wEval,grid, projection)

        print("p =", self._polynomial_degree)
        print("N =", self._num_grid_cells,'  Points per element = ',self._num_eval_points)
        print("maximum error =", max_error)
        
        return
        
    def _dg_plot_mesh(self, zEval,wEval,grid,projection):
    
        exact = self._calculate_exact_solution(grid,self._init_cond,zEval)
        approx = self._calculate_approximate_solution(projection, zEval)
        pointwise_error = exact - approx
        max_error = np.max(np.abs(pointwise_error))

        self._dg_plot_solution_and_approx(grid, exact, approx, self._colors['exact'], self._colors['approx'])
        plt.legend(['Exact', 'Approx'])
        self._dg_plot_semilog_error(grid, np.abs(pointwise_error))
        self._dg_plot_error(grid, exact, approx)

        return max_error
        
    @staticmethod
    def _dg_plot_solution_and_approx(grid, exact, approx, color_exact, color_approx):
        # print(color_exact, color_approx)
        plt.figure(1)
        plt.plot(grid, exact, color_exact)
        plt.plot(grid, approx, color_approx)
        plt.xlabel('x')
        plt.ylabel('$u(x,t)$')
        plt.title('Solution and Approximation')

    @staticmethod
    def _dg_plot_semilog_error(grid, pointwise_error):
        plt.figure(2)
        plt.semilogy(grid[0], pointwise_error[0])
        plt.xlabel('x')
        plt.ylabel('$|u(x,t)-u_h(x,t)|$')
        plt.title('Semilog Error plotted at Evaluation points')

    @staticmethod
    def _dg_plot_error(grid, exact, approx):
        plt.figure(3)
        plt.plot(grid[0], exact[0]-approx[0])
        plt.xlabel('X')
        plt.ylabel('$u(x,t)-u_h(x,t)$')
        plt.title('Errors')
        
    def eval_grid(self):
    
        zEval = np.zeros((self._num_eval_points))
        wEval = np.zeros((self._num_eval_points))
        zEval, wEval = np.polynomial.legendre.leggauss(self._num_eval_points)
        grid = np.zeros((self._num_eval_points*self._num_grid_cells,1))

        for cell in range(self._num_grid_cells):
            eval_points = self._mesh[cell] + 0.5 * self._cell_len * zEval
            eval_points = np.transpose(np.array(eval_points))
           
            if self._num_eval_points == 1:
                grid[cell] = eval_points
            else:
                for j in range(self._num_eval_points):
                    grid[cell*self._num_eval_points+j] = eval_points[j]
            
        return zEval,wEval,grid

    def _calculate_exact_solution(self,grid,init_cond,zEval):
    
        exact = np.zeros((self._num_eval_points*self._num_grid_cells,1))
        for cell in range(len(self._mesh)):
            zpt = np.zeros((self._num_eval_points,1))
            if self._num_eval_points == 1:
                zpt = grid[cell]
            else:
                zpt = grid[cell*self._num_eval_points:(cell+1)*self._num_eval_points-1]
            
            eval_values = init_cond.calculate(zpt)
            if self._num_eval_points == 1:
                exact[cell] = eval_values
            else:
                exact[cell*self._num_eval_points:(cell+1)*self._num_eval_points-1] = eval_values
        
        return exact

    def _calculate_approximate_solution(self, projection, points):
        
        approx = np.zeros((self._num_eval_points*self._num_grid_cells,1))
        basis = self._basis.get_basis_vector()
        

        basis_matrix = [[basis[degree].subs(x, points[point]) for point in range(self._num_eval_points)]
                        for degree in range(self._polynomial_degree+1)]
        
        for cell in range(self._num_grid_cells):
            for pteval in range(self._num_eval_points):
                approx[cell*self._num_eval_points+pteval] = sum(projection[degree][cell] * basis_matrix[degree][pteval] for degree in range(self._polynomial_degree+1))
                  

        return approx
 
########################################################################
    def _reset(self):
        # Set additional necessary instance variables
        #self._right_bound = self._right_bound.astype('float')
        #self._left_bound = self._left_bound.astype('float')
        self._interval_len = self._right_bound-self._left_bound
        #self._interval_len = self._interval_len.astype('float')
        self._cell_len = self._interval_len / self._num_grid_cells
        self._basis = OrthonormalLegendre(self._polynomial_degree)

        # Set mesh with no ghost points on each side
        self._mesh = np.arange(self._left_bound + 0.5*self._cell_len, self._right_bound + (0.5*self._cell_len), self._cell_len)  # +3/2

    def siac_save_plots(self):
        name = self._init_cond.get_name() + '__' +  '__number_of_cells_' \
            + str(self._num_grid_cells) + '__polynomial_degree_' + str(self._polynomial_degree) \
            + '__number_Bsplines_' + str(2*self._RS+1) + '__smoothness_' + str(self._smoothness)

        # Set paths for plot files if not existing already
        if not os.path.exists(self._plot_dir):
            os.makedirs(self._plot_dir)

        if not os.path.exists(self._plot_dir + '/exact_and_filtered'):
            os.makedirs(self._plot_dir + '/exact_and_filtered')

        if not os.path.exists(self._plot_dir + '/semilog_error_filtered'):
            os.makedirs(self._plot_dir + '/semilog_error_filtered')

        if not os.path.exists(self._plot_dir + '/error_filtered'):
            os.makedirs(self._plot_dir + '/error_filtered')

        # Save plots
        plt.figure(4)
        plt.savefig(self._plot_dir + '/exact_and_filtered/' + name + '.pdf')

        plt.figure(5)
        plt.savefig(self._plot_dir + '/semilog_error_filtered/' + name + '.pdf')

        plt.figure(6)
        plt.savefig(self._plot_dir + '/error_filtered/' + name + '.pdf')

        
    def siac_plot_results(self, xEval,filtered,init_cond,zEval):
    
        xEval = np.array(xEval)
        filtered = np.array(filtered)
        max_error = self._siac_plot_mesh(xEval,filtered,init_cond,zEval)

        print("Number of B-Splines = ",2*self._RS+1,"Smoothness = ",self._smoothness)
        print("maximum error =", max_error)
        
    def _siac_plot_mesh(self,grid,filtered,zEval, init_cond):

        exact = self._calculate_exact_solution(grid,self._init_cond,zEval)
        filtered_pointwise_error = np.subtract(exact,filtered)
        
        max_error = np.max(np.abs(filtered_pointwise_error))
        print('Maximum filtered error = ',max_error)

        self._siac_plot_solution_and_filtered(grid, exact, filtered, self._colors['exact'], self._colors['filtered'])
        plt.legend(['Exact', 'Filtered'])
        self._siac_plot_semilog_error(grid, np.abs(filtered_pointwise_error))
        self._siac_plot_error(grid, filtered_pointwise_error)

        return max_error
        
    @staticmethod
    def _siac_plot_solution_and_filtered(grid, exact, filtered, color_exact, color_filtered):
        # print(color_exact, color_approx)
        
        plt.figure(4)
        plt.plot(grid, exact, color_exact,grid, filtered, color_filtered)
#        plt.plot(grid, filtered, color_filtered)
        plt.xlabel('x')
        plt.ylabel('$u(x,t)$')
        plt.title('Solution and Filtered')

    @staticmethod
    def _siac_plot_semilog_error(grid, pointwise_error):
    
        plt.figure(5)
        plt.semilogy(grid, pointwise_error)
        plt.xlabel('x')
        plt.ylabel('$|u(x,t)-u_h^*(x,t)|$')
        plt.title('Semilog Filtered Error plotted at Evaluation points')

    @staticmethod
    def _siac_plot_error(grid, filtered_pointwise_error):

        plt.figure(6)
        plt.plot(grid, filtered_pointwise_error)
        plt.xlabel('X')
        plt.ylabel('$u(x,t)-u_h^*(x,t)$')
        plt.title('Filtered Errors')

