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

# Kullback-Leibler divergence
def kl_div(x,y):
    """
    | :math:`\mbox{kl\_div} :
      \mathbb{R}_{++}^n \\times \mathbb{R}_{++}^n \\to \mathbb{R}, 
      \ \mbox{kl\_div}(x,y) = \sum_{i=1}^n
      ( x_i \mbox{log} (x_i/y_i) - x_i + y_i)`.
    | Convex.

    :param x: number, 
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :param y: number, 
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :return: number or 
             :ref:`tree<tree_obj>`.
    """

    # Check x
    if (np.isscalar(x) or 
        type(x).__name__ in SCALAR_OBJS):
        x = vstack([x])
    elif (type(x) is not cvxpy_matrix and 
          type(x).__name__ not in ARRAY_OBJS):
        raise TypeError('Invalid first argument type')
    elif x.shape[1] != 1:
        raise ValueError('Invalid first argument dimension')
    
    # Check y
    if (np.isscalar(y) or 
        type(y).__name__ in SCALAR_OBJS):
        y = vstack([y])
    elif (type(y) is not cvxpy_matrix and 
          type(y).__name__ not in ARRAY_OBJS):
        raise TypeError('Invalid second argument type')
    elif y.shape[1] != 1:
        raise ValueError('Invalid second argument dimension')

    # Check dimension match
    if x.shape[0] != y.shape[0]:
        raise ValueError('Dimension mismatch')

    # Check positive if numeric
    n = x.shape[0]
    for i in range(0,n,1):
        if np.isscalar(x[i,0]) and x[i,0] <= 0:
            return np.inf
        if np.isscalar(y[i,0]) and y[i,0] <= 0:
            return np.inf

    # Construct objective and constraints
    t = variable()
    u = variable(n,1)
    v = variable(n,1)
    constr = [belongs(vstack((u,v,t)),kl_div_epi),
              geq(u,0),geq(v,0)]
        
    # Construct and return program
    p = program(minimize(t),
                constr,
                [u,v], 
                name='kl_div')
    return p(x,y)
