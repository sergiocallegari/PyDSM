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
from cvxpy.scalars import cvxpy_obj,cvxpy_scalar_var
from cvxpy.utils import vstack,hstack
from cvxpy.interface import belongs
from cvxpy.interface import variable
from cvxpy.interface import greater_equals
from cvxpy.constraints import cvxpy_list
from cvxpy.sets.semidefinite_cone import semidefinite_cone
from functools import reduce

# Class definition
class cvxpy_geo_mean_cone(object):
    """
    | :math:`\{(x,y) \in \mathbb{R}_+^n \\times \mathbb{R} \ |
      \ \\big( \prod_{i=1}^{n}x_i \\big)^{1/n} \ge y \}`.
    """
    
    # Method: __init__
    def __init__(self):
        """
        Class constructor.
        """
        
        self.type = SET
        self.name = 'geo_mean_cone'
        self.expansion_type = PM

    # Method: valid_shape
    def valid_shape(self,shape):

        return shape[0]>=2 and shape[1] == 1

    # Method: __str__
    def __str__(self):

        return self.name

    # Method: _pm_expand
    def _pm_expand(self,constr):
        
        # Get objects
        v = constr.left
        v = vstack([v[i,0] for i in range(0,v.shape[0])])
        n = v.shape[0]-1
        x = vstack([v[0:n,0]])
        y = v[n,0]
        z = variable()

        # One element
        if n==1:
            return cvxpy_list([greater_equals(x,0),
                               greater_equals(x,y)])

        # Get power of 2 size
        m = 0
        while np.log2(n+m) % 1 != 0:
            m += 1

        # Copy elements of x on a list and restrict them
        constr_list = []
        el_list = []
        for i in range(0,n,1):
            el_list+=[x[i,0]]
            if (not np.isscalar(x[i,0]) and
                type(x[i,0]) is not cvxpy_obj):
                constr_list += [greater_equals(x[i,0],0)]

        # Construct expansion
        for i in range(0,m,1):
            el_list += [z]
        while len(el_list) > 2:
            new_list = []
            for i in range(0,int(len(el_list)/2)):
                x1 = el_list[2*i]
                x2 = el_list[2*i+1]
                w = variable()
                constr_list += [belongs(vstack((hstack((x1,w)),
                                                hstack((w,x2)))),
                                        semidefinite_cone)]
                new_list += [w]
            el_list = new_list
        x1 = el_list[0]
        x2 = el_list[1]
        constr_list += [belongs(vstack((hstack((x1,z)),
                                        hstack((z,x2)))),
                                semidefinite_cone)]
        constr_list += [greater_equals(z,0),greater_equals(z,y)]
        return cvxpy_list(constr_list)

    # Method: _construct
    def _construct(self,el,mp,p):
        n = el.shape[0]-1
        u = vstack([el[0:n,0]])
        v = el[n,0]

        # Prepare u
        for i in range(0,n,1):
            if type(u[i,0]) is cvxpy_obj:
                u[i,0] = u[i,0].value
                
        # Prepare v
        if type(v) is cvxpy_obj:
            v = v.value

        # f
        def f(x):
            vals = []
            for i in range(0,n,1):
                if type(u[i,0]) is cvxpy_scalar_var:
                    vals += [x[mp[u[i,0]]]]
                else:
                    vals += [u[i,0]]
            if type(v) is cvxpy_scalar_var:
                y = x[mp[v]]
            else:
                y = v
            return y - reduce(lambda w,z:w*z,vals,1.)**(1./(n*1.))
        
        # grad f
        def grad_f(x):
            g = opt.spmatrix(0.0,[],[],(p,1))
            if type(v) is cvxpy_scalar_var:
                fval = -f(x)+x[mp[v]]
            else:
                fval = -f(x)+v
            for i in range(0,n,1):
                if type(u[i,0]) is cvxpy_scalar_var:
                    xi = x[mp[u[i,0]]]
                    g[mp[u[i,0]]] = -fval/(n*xi*1.)
            if type(v) is cvxpy_scalar_var:
                g[mp[v]] = 1.
            return g

        # hess f
        def hess_f(x):
            h = opt.spmatrix(0.0,[],[],(p,p))
            if type(v) is cvxpy_scalar_var:
                fval = -f(x)+x[mp[v]]
            else:
                fval = -f(x)+v
            for i in range(0,n,1):
                for j in range(0,n,1):
                    if (type(u[i,0]) is cvxpy_scalar_var and 
                        type(u[j,0]) is cvxpy_scalar_var):
                        xi = x[mp[u[i,0]]]
                        xj = x[mp[u[j,0]]]
                        if i == j:
                            w = (n-1.)*fval/(n*n*1.*xi*xi)
                            h[mp[u[i,0]],mp[u[i,0]]] = w
                        else:
                            w = -fval/(1.*n*n*xi*xj)
                            h[mp[u[i,0]],mp[u[j,0]]] = w
                            h[mp[u[j,0]],mp[u[i,0]]] = w
            return h

        # ind f
        def ind_f(x):
            for i in range(0,n,1):
                if type(u[i,0]) is cvxpy_scalar_var:
                    if x[mp[u[i,0]]] <= 0:
                        return False
                else:
                    if u[i,0] <= 0:
                        return False
            return True

        # Return functions
        return f, grad_f, hess_f, ind_f

# Create instance
geo_mean_cone = cvxpy_geo_mean_cone()
