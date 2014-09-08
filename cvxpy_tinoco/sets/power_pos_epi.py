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
# along with this program.  If not, see <http://www.gnu.org/licenses/>. #
#***********************************************************************#

import numpy as np
import cvxopt as opt
from ..defs import *
from ..scalars import cvxpy_obj
from ..scalars import cvxpy_scalar_var

# Class definition
class cvxpy_power_pos_epi(object):
    """
    | :math:`\{ (x,t) \in \mathbb{R}_+ \\times \mathbb{R} \ | 
      \ x^p  \leq t \}` where :math:`p > 1`.
    """

    # Method: __init__
    def __init__(self,p):
        """
        Class constructor.
        """

        self.type = SET
        self.name = 'power_pos_epi'
        self.expansion_type = DIF
        self.p = p

    # Method: valid_shape
    def valid_shape(self,shape):
        
        return shape[1] == 1

    # Method: __str__
    def __str__(self):

        return self.name

    # Method: _construct
    def _construct(self,el,mp,n):

        p = self.p

        # Prepare zg
        if np.isscalar(el[0,0]):
            zg = el[0,0]
        elif type(el[0,0]) is cvxpy_obj:
            zg = el[0,0].value
        else:
            zg = None

        # Prepare tg
        if np.isscalar(el[1,0]):
            tg = el[1,0]
        elif type(el[1,0]) is cvxpy_obj:
            tg = el[1,0].value
        else:
            tg = None
            
        # f
        def f(x):
            z = zg
            t = tg
            if z is None:
                z = x[mp[el[0,0]]]
            if t is None:
                t = x[mp[el[1,0]]]
            return z**p-t

        # grad f
        def grad_f(x):
            z = zg
            g = opt.spmatrix(0.0,[],[],(n,1))
            if z is None:
                z = x[mp[el[0,0]]]
            if type(el[0,0]) is cvxpy_scalar_var:
                g[mp[el[0,0]]] = p*(z**(p-1.))
            if type(el[1,0]) is cvxpy_scalar_var:
                g[mp[el[1,0]]] = -1.
            return g

        # hess f
        def hess_f(x):
            z = zg
            h = opt.spmatrix(0.0,[],[],(n,n))
            if z is None:
                z = x[mp[el[0,0]]]
            if type(el[0,0]) is cvxpy_scalar_var:
                h[mp[el[0,0]],mp[el[0,0]]] = p*(p-1.)*(z**(p-2.))
            return h

        # ind f
        def ind_f(x):
            z = zg
            if z is None:
                z = x[mp[el[0,0]]]
            return z > 0.

        # Return functions
        return f, grad_f, hess_f, ind_f

# Create instance
power_pos_epi = cvxpy_power_pos_epi
