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

# Geometric Mean
def geo_mean(x):
    r"""
    | :math:`\mbox{geo\_mean} :
      \mathbb{R}_+^n \\to \mathbb{R},
      \ \mbox{geo\_mean}(x) = \\big( \prod_{i=1}^n x_i \\big)^{1/n}`.
    | Concave and increasing.

    :param x: number,
              :ref:`scalar object<scalar_ref>` or
              :ref:`multidimensional object<multi_ref>`.
    :return: number or
             :ref:`tree<tree_obj>`.
    """

    # Process input
    if (np.isscalar(x) or
        type(x).__name__ in SCALAR_OBJS):
        x = vstack([x])
    elif (type(x) is not cvxpy_matrix and
          type(x).__name__ not in ARRAY_OBJS):
        raise TypeError('Invalid argument type')

    # x must be a vector
    (m,n) = x.shape
    if n!=1:
        raise ValueError('Invalid argument dimension')

    # Check nonnegativity if cvxpy_matrix
    if (type(x) is cvxpy_matrix and
        not np.all(x >= 0)):
        return -np.inf

    # Construct program
    z = variable(m,1)
    v = variable(m,1)
    t = variable()
    p = program(maximize(t),
                [belongs(vstack((v,t)),geo_mean_cone),
                 greater_equals(z,v), geq(v,0)],
                [z],
                name='geo_mean')
    return p(x)
