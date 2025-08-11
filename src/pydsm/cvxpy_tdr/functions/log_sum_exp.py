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
from ..arrays import cvxpy_array
from ..arrays import cvxpy_matrix

# log_sum_exp
def log_sum_exp(x):
    r"""
    | :math:`\mbox{log\_sum\_exp} :
      \mathbb{R}^n \\to \mathbb{R},
      \ \mbox{log\_sum\_exp}(x) =
      \mbox{log} \\big( \sum_{i = 1}^n e^{x_i} \\big)`.
    | Convex and increasing.

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

    # Must be column
    if x.shape[1] != 1:
        raise ValueError('Invalid argument dimensions')

    # Single element
    if x.shape == (1,1):
        return x[0,0]

    # Construct program
    t = variable()
    z = variable(x.shape[0],1)
    v = variable(x.shape[0],1)
    w = variable(x.shape[0],1)
    constr = [equals(sum(w),1)]
    for i in range(0,x.shape[0],1):
        constr += [belongs(vstack((v[i,0]-t,1.,w[i,0])),
                           exp_cone),
                   less_equals(z[i,0],v[i,0])]
    p = program(minimize(t),
                constr,
                [z],
                name='log_sum_exp')

    # Return
    return p(x)
