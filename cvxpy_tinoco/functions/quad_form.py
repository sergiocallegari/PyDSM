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

# Quadratic form
def quad_form(x,P):
    """
    | :math:`\mbox{quad\_form}:
      \mathbb{R}^n \\times \mathbb{S}_+^n \\to \mathbb{R},
      \ \mbox{quad\_form}(x,P) = x^TPx`.
    | Convex.

    :param x: number,
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :param P: number or
              :ref:`matrix<matrix_obj>`.
    :return: number or 
             :ref:`tree<tree_obj>`.
    """

    # Check first input type
    if (np.isscalar(x) or
        type(x).__name__ in SCALAR_OBJS):
        x = vstack([x])
    elif (type(x) is not cvxpy_matrix and
          type(x).__name__ not in ARRAY_OBJS):
        raise TypeError('Invalid first argument type')
    
    # Check shape
    (m,n) = x.shape
    if n!=1:
        raise ValueError('Invalid first argument dimension')

    # Check second input type
    if np.isscalar(P):
        P *= eye(m)
    elif type(P) is not cvxpy_matrix:
        raise TypeError('Invalid second argument type')

    # Check dimensions
    if P.shape[0] != m:
        raise ValueError('Invalid dimensions')

    # Check symmetry
    if not np.allclose(P,P.T):
        raise ValueError('Invalid second argument')

    # Check positive semidefinite
    if np.min(np.linalg.eig(P)[0]) < -EPSILON:
        raise ValueError('Invalid second argument')

    # Construct and return program
    z = variable(m,1)
    t = variable()
    p = program(minimize(t),
                [belongs(hstack((vstack((eye(m),z.T)),
                                 vstack((z,t)))), 
                         semidefinite_cone)],
                [z], 
                name='quad_form')
    return p(sqrtm(P)*x)
