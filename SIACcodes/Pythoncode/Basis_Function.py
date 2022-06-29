# -*- coding: utf-8 -*-
"""
@author: Laura C. Kühle
modified by:  Jennifer Ryan

"""
import numpy as np
from sympy import Symbol, integrate

x = Symbol('x')
z = Symbol('z')


class Vector(object):
    def __init__(self, polynomial_degree):
        self._polynomial_degree = polynomial_degree
        self._basis = self._build_basis_vector(x)
        self._derivative_basis = self._build_derivative_basis_vector(x)

    def get_basis_vector(self):
        return self._basis

    def get_derivative_basis_vector(self):
        return self._derivative_basis

    def _build_basis_vector(self, eval_point):
        return []

    def _build_derivative_basis_vector(self, eval_point):
        return []

    def get_basis_projections(self):
        pass


class Legendre(Vector):
    def _build_basis_vector(self, eval_point):
        return self._calculate_legendre_vector(eval_point)

    def _build_derivative_basis_vector(self, eval_point):
        return self._calculate_derivative_legendre_vector(eval_point)

    def _calculate_legendre_vector(self, eval_point):
        vector = []
        for degree in range(self._polynomial_degree+1):
            if degree == 0:
                vector.append(1.0 + 0*eval_point)
            else:
                if degree == 1:
                    vector.append(eval_point)
                else:
                    poly = (2.0*degree - 1)/degree * eval_point * vector[-1] - (degree-1)/degree * vector[-2]
                    vector.append(poly)
        return vector

    def _calculate_derivative_legendre_vector(self, eval_point):
        derivative_vector = []
        leg_vector = self._calculate_legendre_vector(eval_point)
        for degree in range(self._polynomial_degree + 1):
            if degree == 0:
                derivative_vector.append(0*eval_point)
            else:
                if degree == 1:
                    derivative_vector.append(1.0 + 0*eval_point)
                else:
                    poly = ((2.0 * degree - 1) / degree) *(leg_vector[degree-1] + (eval_point * derivative_vector[-1])) - ((degree - 1) / degree) * derivative_vector[-2]
                    derivative_vector.append(poly)

        return derivative_vector


class OrthonormalLegendre(Legendre):
    def _build_basis_vector(self, eval_point):
        leg_vector = self._calculate_legendre_vector(eval_point)
        return [leg_vector[degree] * np.sqrt(degree+0.5) for degree in range(self._polynomial_degree+1)]

    def _build_derivative_basis_vector(self, eval_point):
        derivative_leg_vector = self._calculate_derivative_legendre_vector(eval_point)
        return [derivative_leg_vector[degree] * np.sqrt(degree+0.5) for degree in range(self._polynomial_degree+1)]

    

    def get_basis_projections(self):
        basis_projection_left = self._build_basis_matrix(z, 0.5 * (z - 1))
        basis_projection_right = self._build_basis_matrix(z, 0.5 * (z + 1))
        return basis_projection_left, basis_projection_right

    def _build_basis_matrix(self, first_param, second_param):
        matrix = []
        for i in range(self._polynomial_degree + 1):
            row = []
            for j in range(self._polynomial_degree + 1):
                entry = integrate(self._basis[i].subs(x, first_param)
                                  * self._basis[j].subs(x, second_param),
                                  (z, -1, 1))
                row.append(np.float64(entry))
            matrix.append(row)
        return matrix

