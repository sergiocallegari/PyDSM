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
from ..interface import *
from ..arrays import cvxpy_array
from ..arrays import cvxpy_matrix

# max
def max(x):
    """
    | :math:`\mbox{max} :
      \mathbb{R}^{m \\times n} \\to \mathbb{R},
      \ \mbox{max}(X) = \mbox{max} \{ X_{ij} \ | \ 
      i \in [m], \ j \in [n] \}`. 
    | Convex and nondecreasing.

    :param x: number,
              :ref:`scalar object<scalar_ref>` or 
              :ref:`multidimensional object<multi_ref>`.
    :return: number or 
             :ref:`scalar object<scalar_ref>`.
    """

    # Check input
    if (np.isscalar(x) or
        type(x).__name__ in SCALAR_OBJS):
        return x
    elif (type(x) is not cvxpy_matrix and
          type(x).__name__ not in ARRAY_OBJS):
        raise TypeError('Invalid argument type')

    # Single element
    if x.shape == (1,1):
        return x[0,0]

    # Construct program
    t = variable()
    z = variable(x.shape[0],x.shape[1])
    constr = []
    for i in range(0,x.shape[0],1):
        for j in range(0,x.shape[1],1):
            constr += [less_equals(z[i,j],t)]
    p = program(minimize(t),
                constr,
                [z], 
                name='max')

    # Return 
    return p(x)
