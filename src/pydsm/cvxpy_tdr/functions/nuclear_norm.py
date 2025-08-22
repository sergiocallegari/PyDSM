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
from ..sets import *
from ..utils import *
from ..interface import *
from ..arrays import cvxpy_matrix

# Nuclear norm
def nuclear_norm(X):
    r"""
    | :math:`\mbox{nuclear\_norm} :
      \mathbb{R}^{m \\times n} \\to \mathbb{R},
      \ \mbox{nuclear\_norm}(X) = \mbox{Tr} \sqrt{X^TX}`.
    | Convex.

    :param X: number,
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :return: number or
             :ref:`tree<tree_obj>`.
    """

    # Check input
    if (np.isscalar(X) or
        type(X).__name__ in SCALAR_OBJS):
        X = vstack([X])
    elif (type(X) is not cvxpy_matrix and
          type(X).__name__ not in ARRAY_OBJS):
        raise TypeError('Invalid argument type')

    # Construct program
    (m,n) = X.shape
    Y = variable(m,m,'symmetric')
    W = variable(n,n,'symmetric')
    Z = variable(m,n)
    p = program(minimize((1./2.)*(trace(Y)+trace(W))),
                [belongs(vstack((hstack((Y,Z)),
                                 hstack((Z.T,W)))),
                         semidefinite_cone)],
                [Z],
                name='nuclear_norm')
    return p(X)
