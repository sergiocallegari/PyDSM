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
from cvxpy.defs import *
from cvxpy.scalars import cvxpy_obj

# Class definition
class cvxpy_semidefinite_cone(object):
    """
    | :math:`\{X \in \mathbb{R}^{n \\times n} \ |
      \ X = X^T, \ X \succeq 0 \} = \mathbb{S}_+^n`.
    """

    # Method: __init__
    def __init__(self):
        """
        Class constructor.
        """

        self.type = SET
        self.name = 'semidefinite_cone'
        self.expansion_type = SDC
    
    # Method: valid_shape
    def valid_shape(self,shape):
        
        return shape[0] == shape[1]

    # Method: __str__
    def __str__(self):

        return self.name
    
    # Method: _construct
    def _construct(self,el,mp,n):

        m = el.shape[0]
        G = opt.spmatrix(0.0,[],[],(m*m,n))
        h = opt.matrix(0.0,(m*m,1))
        for j in range(0,m,1):
            for i in range(0,m,1):
                if np.isscalar(el[i,j]):
                    h[j*m+i,0] = el[i,j]*1.
                elif type(el[i,j]) is cvxpy_obj:
                    h[j*m+i,0] = el[i,j].value*1.
                else:
                    G[j*m+i,mp[el[i,j]]] = -1.
        return G,h,m

# Create instance
semidefinite_cone = cvxpy_semidefinite_cone()
