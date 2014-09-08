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
from cvxpy.utils import vstack

# Class definition
class cvxpy_kl_div_epi(object):
    """
    | :math:`\{(x,y,t) \in \mathbb{R}_{++}^n \\times \mathbb{R}_{++}^n
      \\times \mathbb{R} \ | 
      \ \sum_{i=1}^n (x_i \mbox{log} (x_i/y_i) - x_i + y_i) \leq t \}`.
    """

    # Method: __init__
    def __init__(self):
        """
        Class constructor.
        """

        self.type = SET
        self.name = 'kl_div_epi'
        self.expansion_type = DIF

    # Method: valid_shape
    def valid_shape(self,shape):
        m,n = shape
        return n == 1 and m >= 3 and m%2 == 1

    # Method: __str__
    def __str__(self):

        return self.name

    # Method: _construct
    def _construct(self,el,mp,p):
        
        n = int((el.shape[0]-1.)/2.)
        u = vstack([el[0:n,0]])
        v = vstack([el[n:2*n,0]])
        t = el[2*n,0]

        # Prepare u and v
        for i in range(0,n):
            if type(u[i,0]) is cvxpy_obj:
                u[i,0] = u[i,0].value
            if type(v[i,0]) is cvxpy_obj:
                v[i,0] = v[i,0].value

        # Prepare t
        if type(t) is cvxpy_obj:
            t = t.value

        # f
        def f(x):
            temp = 0
            for i in range(0,n):
                if type(u[i,0]) is cvxpy_scalar_var:
                    u_i = x[mp[u[i,0]]]
                else:
                    u_i = u[i,0]
                if type(v[i,0]) is cvxpy_scalar_var:
                    v_i = x[mp[v[i,0]]]
                else:
                    v_i = v[i,0]
                temp += u_i*np.log(u_i*1./v_i)-u_i+v_i
            if type(t) is cvxpy_scalar_var:
                return temp - x[mp[t]]
            else:
                return temp - t

        # grad f
        def grad_f(x):
            g = opt.spmatrix(0.0,[],[],(p,1))
            for i in range(0,n):
                if type(u[i,0]) is cvxpy_scalar_var:
                    u_i = x[mp[u[i,0]]]
                else:
                    u_i = u[i,0]
                if type(v[i,0]) is cvxpy_scalar_var:
                    v_i = x[mp[v[i,0]]]
                else:
                    v_i = v[i,0]
                if type(u[i,0]) is cvxpy_scalar_var:
                    g[mp[u[i,0]]] = np.log(u_i*1./v_i)
                if type(v[i,0]) is cvxpy_scalar_var:
                    g[mp[v[i,0]]] = (-u_i*1./v_i)+1.
            if type(t) is cvxpy_scalar_var:
                g[mp[t]] = -1.
            return g

        # hess f
        def hess_f(x):
            h = opt.spmatrix(0.0,[],[],(p,p))
            for i in range(0,n):
                if type(u[i,0]) is cvxpy_scalar_var:
                    u_i = x[mp[u[i,0]]]
                else:
                    u_i = u[i,0]
                if type(v[i,0]) is cvxpy_scalar_var:
                    v_i = x[mp[v[i,0]]]
                else:
                    v_i = v[i,0]
                if type(u[i,0]) is cvxpy_scalar_var:
                    h[mp[u[i,0]],mp[u[i,0]]] = 1./u_i
                if type(v[i,0]) is cvxpy_scalar_var:
                    h[mp[v[i,0]],mp[v[i,0]]] = u_i*1./(v_i**2.)
                if (type(u[i,0]) is cvxpy_scalar_var and
                    type(v[i,0]) is cvxpy_scalar_var):
                    w = -1./v_i
                    h[mp[u[i,0]],mp[v[i,0]]] = w
                    h[mp[v[i,0]],mp[u[i,0]]] = w
            return h

        # ind f
        def ind_f(x):
            for i in range(0,n):
                if type(u[i,0]) is cvxpy_scalar_var:
                    if x[mp[u[i,0]]] <= 0:
                        return False
                else:
                    if u[i,0] <= 0:
                        return False
                if type(v[i,0]) is cvxpy_scalar_var:
                    if x[mp[v[i,0]]] <= 0:
                        return False
                else:
                    if v[i,0] <= 0:
                        return False
            return True

        # Return functions
        return f, grad_f, hess_f, ind_f

# Create instance
kl_div_epi = cvxpy_kl_div_epi()
