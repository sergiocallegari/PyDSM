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
from ..functions.abs import abs
from ..arrays import cvxpy_matrix

# Norm infinity
def norm_inf(X):
    """
    | :math:`\mbox{norm\_inf} :
      \mathbb{R}^{m \\times n} \\to \mathbb{R},
      \ \mbox{norm\_inf}(X) = \mbox{sup} \{ \|Xu\|_{\infty} \ |
      \ \|u\|_{\infty} \leq 1 \}`. 
    | Convex. 

    :param X: number,
              :ref:`scalar object<scalar_ref>` or 
              :ref:`multidimensional object<multi_ref>`.
    :return: number 
             or :ref:`tree<tree_obj>`.
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
    t = variable()
    z = variable(m,n)
    w = vstack([z])                 
    p = program(minimize(t),
                list(map(lambda y: leq(sum(abs(y)),t),
                    [w[i,:] for i in range(0,m)])),
                [z],
                name='norm_inf')
    return p(X)

