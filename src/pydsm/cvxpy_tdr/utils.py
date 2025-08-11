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
from .scalars import cvxpy_obj
from .arrays import cvxpy_array
from .arrays import cvxpy_matrix
from scipy.linalg import sqrtm as sci_sqrtm
from functools import reduce

# Names
__all__ = ["hstack","vstack","sum","randn","rand",
           "seed","eye","zeros","ones","diag","sqrtm", "trace",
           "reshape"]

# Utility function: hstack
def hstack(seq):
    """
    | Stacks objects horizontally.

    :param seq: tuple or list of numbers,
                :ref:`scalar objects<scalar_ref>` or
                :ref:`multidimensional objects<multi_ref>`.
    :return: :ref:`array<array_obj>` or
             :ref:`matrix<matrix_obj>`.
    """

    # Check input
    if (type(seq) is not tuple and
        type(seq) is not list):
        raise TypeError('Invalid argument')

    # Process input
    new_list = []
    numeric = True
    for x in seq:
        if np.isscalar(x):
            new_list += [matrix(x)]
        elif type(x) is cvxpy_obj:
            new_list += [matrix(x.value)]
        elif type(x).__name__ in SCALAR_OBJS:
            numeric = False
            new_x = cvxpy_array(1,1)
            new_x[0,0] = x
            new_list += [new_x]
        elif type(x).__name__ in ARRAY_OBJS:
            numeric = False
            new_list += [x]
        elif type(x) is cvxpy_matrix:
            new_list += [x]
        else:
            raise TypeError('Invalid input')

    # Input is numeric
    if numeric:
        return np.hstack(new_list)

    # Verify dimensions
    m = new_list[0].shape[0]
    for x in new_list:
        if x.shape[0] != m:
            raise ValueError('Invalid dimensions')

    # Allocate new array
    n = int(np.sum([x.shape[1] for x in new_list]))
    new_ar = cvxpy_array(m,n)

    # Fill new array
    k = 0
    for x in new_list:
        for i in range(0,x.shape[0],1):
            for j in range(0,x.shape[1],1):
                new_ar[i,k+j] = x[i,j]
        k = k + x.shape[1]

    # Return new array
    return new_ar

# Utility function: vstack
def vstack(seq):
    """
    | Stacks objects vertically.

    :param seq: tuple or list of numbers,
                :ref:`scalar objects<scalar_ref>` or
                :ref:`multidimensional objects<multi_ref>`.
    :return: :ref:`array<array_obj>` or
             :ref:`matrix<matrix_obj>`.
    """

    # Check input
    if (type(seq) is not tuple and
        type(seq) is not list):
        raise TypeError('Invalid argument')

    # Process input
    new_list = []
    numeric = True
    for x in seq:
        if np.isscalar(x):
            new_list += [matrix(x)]
        elif type(x) is cvxpy_obj:
            new_list += [matrix(x.value)]
        elif type(x).__name__ in SCALAR_OBJS:
            numeric = False
            new_x = cvxpy_array(1,1)
            new_x[0,0] = x
            new_list += [new_x]
        elif type(x).__name__ in ARRAY_OBJS:
            numeric = False
            new_list += [x]
        elif type(x) is cvxpy_matrix:
            new_list += [x]
        else:
            raise TypeError('Invalid Input')

    # Input is numeric
    if numeric:
        return np.vstack(new_list)

    # Verify dimensions
    n = new_list[0].shape[1]
    for x in new_list:
        if x.shape[1] != n:
            raise ValueError('Invalid Dimensions')

    # Allocate new array
    m = int(np.sum([x.shape[0] for x in new_list]))
    new_ar = cvxpy_array(m,n)

    # Fill new array
    k = 0
    for x in new_list:
        for i in range(0,x.shape[0],1):
            for j in range(0,x.shape[1],1):
                new_ar[i+k,j] = x[i,j]
        k = k + x.shape[0]

    # Return new array
    return new_ar

# Utility function: sum
def sum(seq):
    """
    | Sums objects in sequence or entries of array-like object.

    :param seq: number,
                :ref:`scalar object<scalar_ref>`,
                :ref:`multidimensional object<multi_ref>`,
                tuple or list of such objects.
    :return: number,
             :ref:`scalar object<scalar_ref>` or
             :ref:`multidimensional object<multi_ref>`.
    """

    # Number or scalars object
    if (np.isscalar(seq) or
        type(seq).__name__ in SCALAR_OBJS):
        return seq

    # Matrix or object array
    elif (type(seq) is cvxpy_matrix or
          type(seq).__name__ in ARRAY_OBJS):
        return reduce(lambda x,y: x + y,
                      [seq[i,j] for i in range(0,seq.shape[0])
                       for j in range(0,seq.shape[1])],0)

    # List or tuple
    elif type(seq) is tuple or type(seq) is list:
        return reduce(lambda x,y: x + y,seq,0)

    # Other
    else:
        raise TypeError('Invalid argument')

# Utility function: randn
def randn(m,n):
    """
    | See **numpy.random.randn**.

    :param m: rows.
    :param n: columns.
    :return: :ref:`matrix<matrix_obj>`.
    """
    return matrix(np.random.randn(m,n))

# Utility function: rand
def rand(m,n):
    """
    | See **numpy.random.rand**.

    :param m: rows.
    :param n: columns.
    :return: :ref:`matrix<matrix_obj>`.
    """
    return matrix(np.random.rand(m,n))

# Utility function: seed
def seed(s=None):
    """
    | See **numpy.random.seed**.

    :param s: integer, array-like or ``None``.
    :return: ``NoneType``.
    """
    np.random.seed(s)

# Utility function: eye
def eye(n):
    """
    | See **numpy.eye**.

    :param n: rows.
    :return: :ref:`matrix<matrix_obj>`.
    """
    return matrix(np.eye(n))

# Utility function: zeros
def zeros(shape):
    """
    | See **numpy.zeros**.

    :param shape: integer or
                  pair of integers.
    :return: :ref:`matrix<matrix_obj>`.
    """
    return matrix(np.zeros(shape))

# Utility function: ones
def ones(shape):
    """
    | See **numpy.ones**.

    :param shape: integer or
                  pair of integers.
    :return: :ref:`matrix<matrix_obj>`.
    """
    return matrix(np.ones(shape))

# Utility function: diag
def diag(v):
    """
    | Extracts diagonal or constructs a diagonal array.

    :param v: number,
              :ref:`scalar object<scalar_ref>`,
              square or one dimensional
              :ref:`multidimensional object<multi_ref>`.
    :return: number,
             :ref:`scalar object<scalar_ref>` or
             :ref:`multidimensional object<multi_ref>`.
    """

    # Check input
    if (np.isscalar(v) or
        type(v).__name__ in SCALAR_OBJS):
        return v
    elif (type(v) is cvxpy_matrix or
          type(v).__name__ in ARRAY_OBJS):
        if v.shape[1]==1:
            v = v.T
    else:
        raise TypeError('Invalid argument type')

    # Get shape
    (m,n) = v.shape

    # Input is 2D
    if m!=1 and n!=1:

        # Must be square
        if m!=n:
            raise ValueError('Invalid dimensions')

        # cvxpy matrix
        if type(v) is cvxpy_matrix:
            return matrix(np.diag(v)).T

        # Object array
        else:
            new_ar = cvxpy_array(m,1)
            for i in range(0,m,1):
                new_ar[i,0] = v[i,i]
            return new_ar

    # Input is 1D
    else:

        # cvxpy matrix
        if type(v) is cvxpy_matrix:
            return matrix(np.diagflat(v))

        # Object array
        else:
            new_ar = cvxpy_array(n,n)
            for i in range(0,n,1):
                new_ar[i,i] = v[0,i]
            return new_ar

# Utility function: sqrtm
def sqrtm(x):
    """
    | See **scipy.linalg.sqrtm**.

    :param x: positive semidefinite matrix.
    :return: :ref:`matrix<matrix_obj>`.
    """
    return matrix(np.real(sci_sqrtm(x)))

# Utility function: trace
def trace(x):
    """
    | Computes trace.

    :param x: number,
              :ref:`scalar object<scalar_ref>` or
              square :ref:`multidimensional object<multi_ref>`.
    :return: number or
             :ref:`scalar object<scalar_ref>`.
    """

    # Check type
    if (np.isscalar(x) or
        type(x).__name__ in SCALAR_OBJS):
        return x
    elif ((type(x) is not cvxpy_matrix) and
          (type(x).__name__ not in ARRAY_OBJS)):
        raise TypeError('Invalid argument type')

    # Check shape
    if x.shape[0] != x.shape[1]:
        raise ValueError('Invalid dimensions')

    # Get trace
    return sum(diag(x))

# Utility function: reshape
def reshape(v, newshape):
    """
    | Reshapes the array to dimensions newshape (see np.reshape)
    | in FORTRAN (column-major) order.

    :param v: :ref:`array<array_obj>` or
              :ref:`matrix<matrix_obj>`.
    :param newshape: tuple with two integers.
    :return: the reshaped variable or matrix,
             :ref:`array<array_obj>` or
             :ref:`matrix<matrix_obj>`.
    """

    # Check input
    if (type(v) is cvxpy_matrix or
          type(v).__name__ in ARRAY_OBJS):
        pass
    else:
        raise TypeError('Invalid argument type')

    if not (hasattr(newshape, '__iter__') and
            len(newshape) == 2):
        raise TypeError('Invalid argument for newshape')

    # Get shape
    (m,n) = v.shape

    # Get new shape
    (mn,nn) = newshape

    if mn*nn != m*n:
        raise ValueError('Output dimension size does not match input dimension')

    # cvxpy matrix
    if type(v) is cvxpy_matrix:
        return matrix(np.reshape(v, newshape, order='F'))
    # Object array
    else:
        new_ar = cvxpy_array(mn,nn)
        for j in range(0,n,1):
            for i in range(0,m,1):
                k = j*m + i
                new_ar[k % mn,k / mn] = v[i,j]
        return new_ar

# Load modules
from .interface import matrix
