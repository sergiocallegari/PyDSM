#***********************************************************************#
# Copyright (C) 2010-2012 Tomas Tinoco De Rubira                        #
#                                                                       #
# This file is part of CVXPY                                            #     
#                                                                       #
# CVXPY is free software: you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #   
# (at your option) any later version.                                   # 
#                                                                       #
# CVXPY is distributed in the hope that it will be useful,              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with this program.  If not, see <http://www.gnu.org/licenses/>. #
#***********************************************************************#

import numpy as np
from .defs import *
from .arrays import *
from .scalars import *
from .constraints import *

# Names
__all__ = ["variable","variable_reset","parameter",
           "parameter_reset", "matrix","spmatrix","equals","eq",
           "less_equals","leq","greater_equals","geq","program",
           "minimize","maximize","belongs"]

# Interface function: variable
def variable(m=1,n=1,structure=None,name=None):
    """ 
    | Creates an optimization variable.

    :param m: rows.
    :param n: columns.
    :param structure: ``None``,
                      ``'symmetric'``,
                      ``'lower_triangular'`` or
                      ``'upper_triangular'``.
    :param name: string or 
                 ``None``.
    :return: :ref:`scalar variable<scalar_var_obj>` or
             :ref:`multidimensional variable<multi_var_obj>`.
    """

    if (m,n) == (1,1):
        return cvxpy_scalar_var(name)
    else:
        return cvxpy_var(m,n,structure,name)

# Interface function: variable_reset
def variable_reset():
    """
    | Resets variable counter.
    """

    cvxpy_scalar_var.i = 0

# Interface function: parameter
def parameter(m=1,n=1,attribute=None,name=None):
    """ 
    | Creates a parameter.

    :param m: rows.
    :param n: columns.
    :param attribute: ``None``,
                      ``'nonnegative'`` or
                      ``'nonpositive'``.
    :param name: string or 
                 ``None``.
    :return: :ref:`scalar parameter<scalar_param_obj>` or
             :ref:`multidimensional parameter<multi_param_obj>`.
    """

    if (m,n) == (1,1):
        return cvxpy_scalar_param(attribute,name)
    else:
        return cvxpy_param(m,n,attribute,name)

# Interface function: parameter_reset
def parameter_reset():
    """
    | Resets parameter counter.
    """

    cvxpy_scalar_param.i = 0

# Interface function: matrix
def matrix(data):
    """ 
    | Creates a matrix. See **numpy.matrix**.
    
    :param data: array-like or string.
    :return: :ref:`matrix<matrix_obj>`.
    """

    # Check if data is numpy matrix
    if type(data) is np.matrix:
        data = np.array(data)

    # Return matrix
    return cvxpy_matrix(data,np.float64)

# Interface function: spmatrix
def spmatrix(data):
    """ 
    | Creates a sparse matrix. See **scipy.sparse.lil_matrix**.
    
    :param data: Dense matrix, sparse matrix or shape tuple.
    :return: :ref:`sparse matrix<spmatrix_obj>`.
    """

    # Return sparse matrix
    return cvxpy_spmatrix(data,np.float64)

# Interface function: equals
def equals(l,r):
    """
    | Forms constraint :math:`l = r`.

    :param l: number, 
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :param r: number, 
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :return: :ref:`constraint<constr_obj>` or
             :ref:`list of constraints<constr_list_obj>`.
    """
    return compare(l,EQUALS,r)

# Interface function: eq
eq = equals

# Interface function: less_equals
def less_equals(l,r):
    """
    | Forms constraint :math:`l \leq r`.

    :param l: number, 
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :param r: number, 
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :return: :ref:`constraint<constr_obj>` or
             :ref:`list of constraints<constr_list_obj>`.
    """
    return compare(l,LESS_EQUALS,r)

# Interface function: leq
leq = less_equals

# Interface function: greater_equals
def greater_equals(l,r):
    """
    | Forms constraint :math:`l \geq r`.

    :param l: number, 
              :ref:`scalar object<scalar_ref>` or 
              :ref:`multidimensional object<multi_ref>`.
    :param r: number, 
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :return: :ref:`constraint<constr_obj>` or
             :ref:`list of constraints<constr_list_obj>`.
    """
    return compare(l,GREATER_EQUALS,r)

# Interface function: geq
geq = greater_equals

# Interface function: program
def program(pair,constraints=[],formals=[],options=None,name='prog'):
    """
    | Creates an optimization program.

    :param pair: output of function 
                 :func:`minimize() <cvxpy.interface.minimize>` or
                 :func:`maximize() <cvxpy.interface.maximize>`.
    :param constraints: list of :ref:`constraints<constr_obj>`.
    :param formals: list of 
                    :ref:`scalar variables<scalar_var_obj>`,
		    :ref:`scalar parameters<scalar_param_obj>`,
                    :ref:`multidimensional variables<multi_var_obj>` or
		    :ref:`multidimensional parameters<multi_param_obj>`.
    :param options: dictionary or
                    ``None``.
    :param name: string.
    :return: :ref:`program<program_obj>`.
    """

    # Get action and objective function
    action = pair[0]
    objective = pair[1]

    # Create and return program
    return cvxpy_program(action,objective,constraints,formals,
                         options,name)

# Interface function: minimize
def minimize(objective):
    """
    | Forms pair (action, objective function).
    
    :param objective: number or 
                      :ref:`scalar object<scalar_ref>`.
    :return: tuple.
    """
    return MINIMIZE,objective

# Interface function: maximize
def maximize(objective):
    """
    | Forms pair (action, objective function).
    
    :param objective: number or 
                      :ref:`scalar object<scalar_ref>`.
    :return: tuple.
    """
    return MAXIMIZE,objective

# Method: belongs
def belongs(x,S):
    """
    | Forms constraint :math:`x \in S`.
    
    :param x: :ref:`multidimensional object<multi_ref>`.
    :param S: set instance.
    :rtype: :ref:`constraint<constr_obj>`.
    """

    # Check first argument
    if ((type(x) is not cvxpy_matrix) and
        (type(x).__name__ not in ARRAY_OBJS)):
        raise TypeError('Invalid first argument')

    # Convert cvxpy_matrix to array
    if type(x) is cvxpy_matrix:
        new_ar = cvxpy_array(x.shape[0],x.shape[1])
        for i in range(0,x.shape[0],1):
            for j in range(0,x.shape[1],1):
                new_ar[i,j] = x[i,j]
        x = new_ar

    # Check second argument
    try:
        if S.type != SET:
            raise TypeError('Invalid second argument')
    except AttributeError:
        raise TypeError('Invalid second argument')

    # Verify dimensions
    if not S.valid_shape(x.shape):
        raise ValueError('Invalid dimensions')

    # Construct and return constraint
    return cvxpy_constr(x,BELONGS,S)

# Load modules
from .programs import *
