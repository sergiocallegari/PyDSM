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
from ..sets import *
from ..utils import *
from ..interface import *
from ..arrays import cvxpy_matrix

# Lambda max
def lambda_max(x):
    """
    | :math:`\mbox{lambda\_max} :
      \mathbb{S}^n \\to \mathbb{R},
      \ \mbox{lambda\_max}(X) = \mbox{sup}\{ y^TXy \ | \ \|y\|_2 = 1 \}`.
    | Convex.
    
    :param x: number,
              :ref:`scalar object<scalar_ref>` or 
              :ref:`multidimensional object<multi_ref>`.
    :return: number or 
             :ref:`tree<tree_obj>`.
    """
    
    # Check type
    if (np.isscalar(x) or
        type(x).__name__ in SCALAR_OBJS):
        return x
    elif (type(x) is not cvxpy_matrix and
          type(x).__name__ not in ARRAY_OBJS):
        raise TypeError('Invalid argument type')

    # Check shape
    (m,n) = x.shape
    if m!=n:
        raise ValueError('Invalid argument dimensions')
    
    # Construct objective and constraints
    z = variable(m,m)
    t = variable()
    obj = t
    constr = [belongs(t*eye(n)-z,semidefinite_cone)]

    # Add symmetry constraint if input is not numeric
    if type(x).__name__ in ARRAY_OBJS:
        constr += [equals(z,z.T)]

    # Check symmetry if input is numeric
    elif not np.allclose(x,x.T):
        return np.inf
     
    # Construct program
    p = program(minimize(obj),
                constr,
                [z],
                name='lambda_max')

    # Return program
    return p(x)
