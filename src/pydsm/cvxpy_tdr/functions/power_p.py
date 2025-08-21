#***********************************************************************#
# Copyright (C) 2010-2013 Tomas Tinoco De Rubira                        #
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>. #
#***********************************************************************#

import numpy as np
from ..defs import *
from ..utils import *
from ..sets import *
from ..interface import *
from ..arrays import cvxpy_matrix
from ..arrays import cvxpy_array

# power_p
def power_p(x,p):
    r"""
    | :math:`\mbox{power\_p} :
      \mathbb{R}^{m \\times n} \\times [1,\infty) \\to \mathbb{R}^{m \\times n},
      \ \mbox{power\_p}(X,p)_{ij} = \\left\{
      \\begin{array}{ll}
      X_{ij}^p, & X_{ij} \geq 0 \\\\
      +\infty,  & X_{ij} < 0.
      \\end{array} \\right.`
    | Convex.

    :param x: number,
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :param p: number, greater than or equal to one.
    :return: number,
             :ref:`tree<tree_obj>`,
             :ref:`matrix<matrix_obj>` or
             :ref:`array<array_obj>`.
    """

    # Prepare input
    if (np.isscalar(x) or
        type(x).__name__ in SCALAR_OBJS):
        arg = vstack([x])
    elif (type(x) is cvxpy_matrix or
          type(x).__name__ in ARRAY_OBJS):
        arg = x
    else:
        raise TypeError('Invalid argument')

    # Check p
    if not (np.isscalar(p)):
        raise TypeError('Invalid argument')
    elif p < 1.:
        raise ValueError('Invalid argument value')

    # Prepare output
    if type(arg) is cvxpy_matrix:
        output = zeros(arg.shape)
    else:
        output = cvxpy_array(arg.shape[0],arg.shape[1])

    # Construct program
    for i in range(0,arg.shape[0],1):
        for j in range(0,arg.shape[1],1):
            t = variable()
            z = variable()
            if p > 1.:
                prog = program(minimize(t),
                               [belongs(vstack((z,t)),
                                        power_pos_epi(p)),
                                greater_equals(z,0)],
                               [z],
                               name='power_p')
            else:
                prog = program(minimize(t),
                               [less_equals(z,t),
                                greater_equals(z,0)],
                               [z],
                               name='power_p')
            if np.isscalar(arg[i,j]) and arg[i,j] < 0.:
                output[i,j] = np.inf
            else:
                output[i,j] = prog(arg[i,j])

    # Return output
    if output.shape == (1,1):
        return output[0,0]
    else:
        return output
