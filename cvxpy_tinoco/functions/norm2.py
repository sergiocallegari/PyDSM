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
from ..functions.lambda_max import lambda_max

# Norm 2
def norm2(X):
    """
    | :math:`\mbox{norm2} :
      \mathbb{R}^{m \\times n} \\to \mathbb{R},
      \ \mbox{norm2}(X) = \mbox{sup} \{ \|Xu\|_2 \ | \ \|u\|_2 \leq 1 \}`.
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
    elif (type(X) is cvxpy_matrix or
          type(X).__name__ in ARRAY_OBJS):
        X = X.T if X.shape[0] == 1 else X
    else:
        raise TypeError('Invalid argument type')

    # Get shape
    (m,n) = X.shape
    
    # Vector
    if 1 in (m,n):
        t = variable()
        z = variable(m,1)
        p = program(minimize(t),
                    [belongs(vstack((z,t)),second_order_cone)],
                    [z],
                    name='norm2')
    # Matrix
    else:
        t = variable()
        z = variable(m,n)
        wm = zeros((m,m))
        wn = zeros((n,n))
        h = vstack([hstack([wm,z]),
                    hstack([z.T,wn])])
        p = program(minimize(t),
                    [leq(lambda_max(h),t)],
                    [z],
                    name='norm2')
    return p(X)
