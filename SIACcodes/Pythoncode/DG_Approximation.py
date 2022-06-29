# -*- coding: utf-8 -*-
"""
@author: Laura C. Kühle

Modified for SIAC by:  Jennifer Ryan

"""
import numpy as np
from sympy import Symbol
<<<<<<< HEAD

import math
#import matplotlib.pyplot as plt
=======
import os
import math
import matplotlib.pyplot as plt
>>>>>>> 00ce8dbc94665b55f9b3feea2985aa98ad0b479e


import Initial_Condition
#import Quadrature
#import OutputInfo
from Basis_Function import OrthonormalLegendre

x = Symbol('x')
<<<<<<< HEAD
z = Symbol('z')
=======

>>>>>>> 00ce8dbc94665b55f9b3feea2985aa98ad0b479e

class DGScheme(object):
    """
    Do documentation here.

    Here come some parameter.
    """

    def __init__(self, **kwargs):
        # Unpack keyword arguments
        self._polynomial_degree = kwargs.pop('polynomial_degree', 2)
        self._num_grid_cells = kwargs.pop('num_grid_cells', 64)
        self._verbose = kwargs.pop('verbose', False)
        self._left_bound = kwargs.pop('left_bound', -1)
        self._right_bound = kwargs.pop('right_bound', 1)
        self._num_quad_points = self._polynomial_degree + 1
<<<<<<< HEAD
#        self._num_eval_points = kwargs.pop('num_eval_points', 1)
=======
        self._num_eval_points = kwargs.pop('num_eval_points', 6)
>>>>>>> 00ce8dbc94665b55f9b3feea2985aa98ad0b479e
        self._plot_dir = kwargs.pop('plot_dir', 'fig')
        self._colors = kwargs.pop('colors', {})
        
        self._init_cond = kwargs.pop('init_cond', 'Sine')
        self._init_config = kwargs.pop('init_config', {})
        #self._quadrature = kwargs.pop('quadrature', 'Gauss')
        #self._quadrature_config = kwargs.pop('quadrature_config', {})
        #self._output = kwargs.pop('output','OutputInfo')
        #self._output_config = kwargs.pop('output_config',{})
        
        
        # Throw an error if there are extra keyword arguments
        if len(kwargs) > 0:
            extra = ', '.join('"%s"' % k for k in list(kwargs.keys()))
            raise ValueError('Unrecognized arguments: %s' % extra)

        # Make sure all classes actually exist
        if not hasattr(Initial_Condition, self._init_cond):
            raise ValueError('Invalid initial condition: "%s"' % self._init_cond)
        #if not hasattr(Quadrature, self._quadrature):
        #    raise ValueError('Invalid quadrature: "%s"' % self._quadrature)
        #if not hasattr(OutputInfo, self._output):
        #    raise ValueError('Invalid output: "%s"' % self._output)
        
        self._check_colors()
        self._reset()

        # Replace the string names with the actual class instances
        # (and add the instance variables for the quadrature)
        self._init_cond = getattr(Initial_Condition, self._init_cond)(self._left_bound, self._right_bound, self._init_config)
        #self._quadrature = getattr(Quadrature, self._quadrature)(self._quadrature_config)
        #self._output = getattr(OutputInfo, self._output)(self._output_config)
    
    def _check_colors(self):
        self._colors['exact'] = self._colors.get('exact', 'k-')
        self._colors['approx'] = self._colors.get('approx', 'r--')

    def approximate(self):
        """
        Do documentation here.

        Here come some parameter.
        """

        projection = self._do_initial_projection(self._init_cond)
        # Plot exact/approximate results, errors, shock tubes and any detector-dependant plots
        #self._output.plot_results(projection, self._init_cond.__class__.__name__)
        #self._output.plot_results(projection)
        
        
        #approx = self._calculate_approximate_solution(projection[:, 1:-1], zEval, self._polynomial_degree)
        
        #if self._verbose:
        #    plt.show()
            
        return projection

<<<<<<< HEAD
=======
    def _calculate_approximate_solution(self, projection, points):
        num_points = len(points)
        basis = self._basis.get_basis_vector()

        basis_matrix = [[basis[degree].subs(x, points[point]) for point in range(num_points)]
                        for degree in range(self._polynomial_degree+1)]

        approx = [[sum(projection[degree][cell] * basis_matrix[degree][point]
                       for degree in range(self._polynomial_degree+1))
                   for point in range(num_points)]
                  for cell in range(len(projection[0]))]

        return np.reshape(np.array(approx), (1, len(approx) * num_points))
>>>>>>> 00ce8dbc94665b55f9b3feea2985aa98ad0b479e

    def _reset(self):
        # Set additional necessary instance variables
        #self._right_bound = self._right_bound.astype('float')
        #self._left_bound = self._left_bound.astype('float')
        self._interval_len = self._right_bound-self._left_bound
        #self._interval_len = self._interval_len.astype('float')
        self._cell_len = self._interval_len / self._num_grid_cells
        self._basis = OrthonormalLegendre(self._polynomial_degree)

        # Set mesh with one ghost point on each side
        self._mesh = np.arange(self._left_bound - (3/2*self._cell_len), self._right_bound + (5/2*self._cell_len), self._cell_len)  # +3/2

        # Set inverse mass matrix
        mass_matrix = []
        for i in range(self._polynomial_degree+1):
            new_row = []
            for j in range(self._polynomial_degree+1):
                new_entry = 0.0
                if i == j:
                    new_entry = 1.0
                new_row.append(new_entry)
            mass_matrix.append(new_row)
        self._inv_mass = np.array(mass_matrix)

    def _do_initial_projection(self, initial_condition, adjustment=0):
        # Initialize matrix and set first entry to accommodate for ghost cell
<<<<<<< HEAD
        output_matrix = np.zeros((self._polynomial_degree+1,self._num_grid_cells))
=======
        output_matrix = [0]
>>>>>>> 00ce8dbc94665b55f9b3feea2985aa98ad0b479e
        basis_vector = self._basis.get_basis_vector()
        
        #THERE IS A DIFFERENCE BETWEEN THE QUADRATURE POINTS TO CALCULATE THE MODES AND THE QUADRATURE POINTS TO EVALUATE THE APPROXIMATION.
        zquad = np.zeros((self._num_quad_points))
        wquad = np.zeros((self._num_quad_points))
        # Gauss-Legendre (default interval is [-1, 1])
        zquad, wquad = np.polynomial.legendre.leggauss(self._num_quad_points)

        for cell in range(self._num_grid_cells):
<<<<<<< HEAD
            # eval_point are the cell centers
            eval_point = self._left_bound + (cell+0.5)*self._cell_len

            for degree in range(self._polynomial_degree + 1):
                new_entry = sum(initial_condition.calculate(
                        eval_point + 0.5 * self._cell_len * zquad[point] - adjustment)
                    * basis_vector[degree].subs(x, zquad[point])
                    * wquad[point] for point in range(self._num_quad_points))
                output_matrix[degree,cell] = new_entry

        return output_matrix

    

        
    
=======
            new_row = []
            eval_point = self._left_bound + (cell+0.5)*self._cell_len

            for degree in range(self._polynomial_degree + 1):
                #new_entry = sum(
                    #initial_condition.calculate(
                    #    eval_point + self._cell_len/2 * self._quadrature.get_quad_points()[point] - adjustment)
                    #* basis_vector[degree].subs(x, self._quadrature.get_quad_points()[point])
                    #* self._quadrature.get_quad_weights()[point]
                new_entry = sum(initial_condition.calculate(
                        eval_point + self._cell_len/2 * zquad[point] - adjustment)
                    * basis_vector[degree].subs(x, zquad[point])
                    * wquad[point] for point in range(self._num_quad_points))
                new_row.append(np.float64(new_entry))

            new_row = np.array(new_row).T
            output_matrix.append(self._inv_mass @ new_row)

        # Set ghost cells to respective value
        output_matrix[0] = output_matrix[self._num_grid_cells]
        output_matrix.append(output_matrix[1])

        # print(np.array(output_matrix).shape)
        return np.transpose(np.array(output_matrix))

    def save_plots(self):
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

        
    def plot_results(self, projection):
        max_error = self._plot_mesh(projection)

        print("p =", self._polynomial_degree)
        print("N =", self._num_grid_cells)
        print("maximum error =", max_error)
        
    def _plot_mesh(self, projection):
    
        zEval = np.zeros((self._num_eval_points))
        wEval = np.zeros((self._num_eval_points))
        zEval, wEval = np.polynomial.legendre.leggauss(self._num_eval_points)
        
        grid, exact = self._calculate_exact_solution(self._mesh[2:-2], self._cell_len, self._init_cond,zEval)
        approx = self._calculate_approximate_solution(projection[:, 1:-1], zEval)

        pointwise_error = np.abs(exact-approx)
        max_error = np.max(pointwise_error)

        self._plot_solution_and_approx(grid, exact, approx, self._colors['exact'], self._colors['approx'])
        plt.legend(['Exact', 'Approx'])
        self._plot_semilog_error(grid, pointwise_error)
        self._plot_error(grid, exact, approx)

        return max_error
        
    @staticmethod
    def _plot_solution_and_approx(grid, exact, approx, color_exact, color_approx):
        # print(color_exact, color_approx)
        plt.figure(1)
        plt.plot(grid[0], exact[0], color_exact)
        plt.plot(grid[0], approx[0], color_approx)
        plt.xlabel('x')
        plt.ylabel('u(x,t)')
        plt.title('Solution and Approximation')

    @staticmethod
    def _plot_semilog_error(grid, pointwise_error):
        plt.figure(2)
        plt.semilogy(grid[0], pointwise_error[0])
        plt.xlabel('x')
        plt.ylabel('|u(x,t)-uh(x,t)|')
        plt.title('Semilog Error plotted at Evaluation points')

    @staticmethod
    def _plot_error(grid, exact, approx):
        plt.figure(3)
        plt.plot(grid[0], exact[0]-approx[0])
        plt.xlabel('X')
        plt.ylabel('u(x,t)-uh(x,t)')
        plt.title('Errors')

    def _calculate_exact_solution(self, mesh, cell_len, initial_condition, zEval):
        grid = []
        exact = []

        for cell in range(len(mesh)):
            eval_points = mesh[cell] + cell_len/2 * zEval

            eval_values = []
            for point in range(len(eval_points)):
                new_entry = self._init_cond.calculate(eval_points[point])
                eval_values.append(new_entry)

            grid.append(eval_points)
            exact.append(eval_values)

        exact = np.reshape(np.array(exact), (1, len(exact) * len(exact[0])))
        grid = np.reshape(np.array(grid), (1, len(grid) * len(grid[0])))

        return grid, exact
>>>>>>> 00ce8dbc94665b55f9b3feea2985aa98ad0b479e
