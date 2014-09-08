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
from cvxpy.defs import *
from cvxpy.sets import *
from cvxpy.utils import *
from cvxpy.interface import *
from cvxpy.arrays import cvxpy_matrix
from cvxpy.functions.geo_mean import geo_mean

# det_rootn
def det_rootn(x):
    """
    | :math:`\mbox{det\_rootn} : 
      \mathbb{S}_+^n \\to \mathbb{R},
      \ \mbox{det\_rootn}(X) = (\mbox{det}(X))^{1/n}`.
    | Concave.

    :param x: number,
              :ref:`scalar object<scalar_ref>` or 
              :ref:`multidimensional object<multi_ref>`.
    :return: number or 
             :ref:`tree<tree_obj>`.
    """

    # Check type
    if (np.isscalar(x) or 
        type(x).__name__ in SCALAR_OBJS):
        x = vstack([x])
    elif (type(x) is not cvxpy_matrix and
          type(x).__name__ not in ARRAY_OBJS):
        raise TypeError('Invalid argument type')

    # Check shape
    (m,n) = x.shape
    if m!=n:
        raise ValueError('Invalid argument dimensions')

    # Construct objective and constraints
    v = variable(m,m)
    Z = variable(m,m,'lower_triangular')
    D = diag(Z)

    constr = [belongs(vstack((hstack((diag(D),Z.T)),
                              hstack((Z,v)))),
                      semidefinite_cone),geq(D,0)]
    
    # Add constraints if argument is not numeric
    if type(x).__name__ in ARRAY_OBJS:
        constr += [belongs(vstack([v]),semidefinite_cone),
                   equals(v,v.T)]
        
    # Check properties if argument is nunmeric
    elif (not np.allclose(x,x.T) or
          np.min(np.linalg.eig(x)[0]) < 0):
        return -np.inf

    # Construct program
    obj = geo_mean(D)
    p = program(maximize(obj),
                constr,
                [v],
                name='det_rootn')

    # Return program
    return p(x)
         
