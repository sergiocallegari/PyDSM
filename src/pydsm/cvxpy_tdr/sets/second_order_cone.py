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
import cvxopt as opt
from ..defs import *
from ..scalars import cvxpy_obj

# Class definition
class cvxpy_second_order_cone(object):
    r"""
    | :math:`\{(x,y) \in \mathbb{R}^n \\times \mathbb{R} \ |
      \ \|x\|_2 \le y \}`.
    """

    # Method: __init__
    def __init__(self):
        """
        Class constructor.
        """

        self.type = SET
        self.name = 'second_order_cone'
        self.expansion_type = SOC

    # Method: valid_shape
    def valid_shape(self,shape):

        return shape[0] > 1 == shape[1]

    # Method: __str__
    def __str__(self):

        return self.name

    # Method: _construct
    def _construct(self,el,mp,n):

        m = el.shape[0]-1
        G = opt.spmatrix(0.0,[],[],(m+1,n))
        h = opt.matrix(0.0,(m+1,1))

        # y
        y = el[m,0]
        if np.isscalar(y):
            h[0,0] = y*1.
        elif type(y) is cvxpy_obj:
            h[0,0] = y.value*1.
        else:
            G[0,mp[y]]=-1.

        # x
        for i in range(0,m,1):
            x = el[i,0]
            if np.isscalar(x):
                h[i+1,0] = x*1.
            elif type(x) is cvxpy_obj:
                h[i+1,0] = x.value*1.
            else:
                G[i+1,mp[x]] = -1.

        # Return G,h
        return G,h,m+1

# Create instance
second_order_cone = cvxpy_second_order_cone()
