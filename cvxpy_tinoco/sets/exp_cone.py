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
from cvxpy.scalars import cvxpy_scalar_var

# Class definition
class cvxpy_exp_cone(object):
    """
    | :math:`\{(x,y,z) \in 
      \mathbb{R} \\times \mathbb{R}_{++} \\times \mathbb{R}
      \ | \ ye^{x/y} \le z \}`.
    """

    # Method: __init__
    def __init__(self):
        """
        Class constructor.
        """

        self.type = SET
        self.name = 'exp_cone'
        self.expansion_type = DIF

    # Method: valid_shape
    def valid_shape(self,shape):
        
        return shape == (3,1)

    # Method: __str__
    def __str__(self):

        return self.name

    # Method: _construct
    def _construct(self,el,mp,n):
        
        # Prepare vg
        if np.isscalar(el[0,0]):
            vg = el[0,0]
        elif type(el[0,0]) is cvxpy_obj:
            vg = el[0,0].value
        else:
            vg = None

        # Prepare yg
        if np.isscalar(el[1,0]):
            yg = el[1,0]
        elif type(el[1,0]) is cvxpy_obj:
            yg = el[1,0].value
        else:
            yg = None
            
        # Prepare zg
        if np.isscalar(el[2,0]):
            zg = el[2,0]
        elif type(el[2,0]) is cvxpy_obj:
            zg = el[2,0].value
        else:
            zg = None

        # f
        def f(x):
            v = vg
            y = yg
            z = zg
            if v is None:
                v = x[mp[el[0,0]]]
            if y is None:
                y = x[mp[el[1,0]]]
            if z is None:
                z = x[mp[el[2,0]]]
            return y*np.exp((v*1.)/(y*1.))-z

        # grad f
        def grad_f(x):
            v = vg
            y = yg
            g = opt.spmatrix(0.0,[],[],(n,1))
            if v is None:
                v = x[mp[el[0,0]]]
            if y is None:
                y = x[mp[el[1,0]]]
            if type(el[0,0]) is cvxpy_scalar_var:
                g[mp[el[0,0]]] = np.exp((v*1.)/(y*1.))
            if type(el[1,0]) is cvxpy_scalar_var:
                t = (v*1.)/(y*1.)
                g[mp[el[1,0]]] = np.exp(t)*(1-t)
            if type(el[2,0]) is cvxpy_scalar_var:
                g[mp[el[2,0]]] = -1.
            return g

        # hess f
        def hess_f(x):
            v = vg
            y = yg
            h = opt.spmatrix(0.0,[],[],(n,n))
            if v is None:
                v = x[mp[el[0,0]]]
            if y is None:
                y = x[mp[el[1,0]]]
            if type(el[0,0]) is cvxpy_scalar_var:
                h[mp[el[0,0]],mp[el[0,0]]] = np.exp((v*1.)/(y*1.))/(y*1.)
            if (type(el[0,0]) is cvxpy_scalar_var and
                type(el[1,0]) is cvxpy_scalar_var):
                t = (v*1.)/(y*1.)
                w = -np.exp(t)*t/(y*1.)
                h[mp[el[0,0]],mp[el[1,0]]] = w
                h[mp[el[1,0]],mp[el[0,0]]] = w
            if type(el[1,0]) is cvxpy_scalar_var:
                t = (v*1.)/(y*1.)
                h[mp[el[1,0]],mp[el[1,0]]] = np.exp(t)*(t**2.)/(y*1.)
            return h

        # ind f
        def ind_f(x):
            y = yg
            if y is None:
                return x[mp[el[1,0]]] > 0
            else:
                return y > 0

        # Return functions
        return f, grad_f, hess_f, ind_f

# Create instance
exp_cone = cvxpy_exp_cone()
