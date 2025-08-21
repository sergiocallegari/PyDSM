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
# along with this program.  If not, see <https://www.gnu.org/licenses/>. #
#***********************************************************************#

import numpy as np
from ..defs import *
from ..utils import *
from ..sets import *
from ..interface import *
from ..arrays import cvxpy_matrix
from ..arrays import cvxpy_array

# square
def square(x):
    r"""
    | :math:`\mbox{square} :
      \mathbb{R}^{m \\times n} \\to \mathbb{R}^{m \\times n},
      \ \mbox{square}(X)_{ij} = X_{ij}^2`.
    | Convex.

    :param x: number,
             :ref:`scalar object<scalar_ref>` or
             :ref:`multidimensional object<multi_ref>`.
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
            p = program(minimize(t),
                        [belongs(vstack((hstack((1.,z)),
                                         hstack((z,t)))),
                                 semidefinite_cone)],
                        [z],
                        name='square')
            output[i,j] = p(arg[i,j])

    # Return output
    if output.shape == (1,1):
        return output[0,0]
    else:
        return output
