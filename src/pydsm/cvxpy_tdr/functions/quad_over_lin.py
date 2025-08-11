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
from ..defs import *
from ..utils import *
from ..sets import *
from ..interface import *
from ..arrays import cvxpy_matrix

# Quadratic over linear
def quad_over_lin(x,y):
    r"""
    | :math:`\mbox{quad\_over\_lin} :
      \mathbb{R}^n \\times \mathbb{R}_{++} \\to \mathbb{R},
      \ \mbox{quad\_over\_lin}(x,y) = x^Tx/y`.
    | Convex in :math:`(x,y)` and decreasing in :math:`y`.

    :param x: number,
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :param y: number or
              :ref:`scalar object<scalar_ref>`.
    :return: number or
             :ref:`tree<tree_obj>`.
    """

    # Check x
    if (np.isscalar(x) or
        type(x).__name__ in SCALAR_OBJS):
        x = vstack([x])
    elif (type(x) is not cvxpy_matrix and
          type(x).__name__ not in ARRAY_OBJS):
        raise TypeError('Invalid first argument')
    elif x.shape[1] != 1:
        raise TypeError('Invalid first argument')

    # Check y
    if (not np.isscalar(y) and
        type(y).__name__ not in SCALAR_OBJS):
        raise TypeError('Invalid second argument')

    # Check that y>0 if y is numeric
    if np.isscalar(y) and y <= 0:
        return np.inf

    # Construct objective and constraints
    n = x.shape[0]
    z = variable(n,1)
    w = variable()
    v = variable()
    t = variable()
    A = vstack((hstack((diag(v*ones((1,n))),z)),
                hstack((z.T,t))))
    constr = [belongs(A,semidefinite_cone),
              greater_equals(w,v)]

    # Construct and return program
    p = program(minimize(t),
                constr,
                [z,w],
                name='quad_over_lin')
    return p(x,y)
