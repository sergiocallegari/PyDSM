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
from ..scalars import cvxpy_scalar_var
from scipy.stats import norm

# Class definition
class cvxpy_log_norm_cdf_hypo(object):
    """
    | :math:`\{ (x,t) \in \mathbb{R} \\times \mathbb{R} \ | 
      \ \mbox{log} ( \int_{-\infty}^x \\frac{1}{\sqrt{2\pi}} 
      e^{-\\frac{u^2}{2}} \mathrm{d}u ) \geq t \}`.
    """

    # Method: __init__
    def __init__(self):
        """
        Class constructor.
        """

        self.type = SET
        self.name = 'log_norm_cdf_hypo'
        self.expansion_type = DIF

    # Method: valid_shape
    def valid_shape(self,shape):
        
        return shape == (2,1)

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
            
        # f
        def f(x):
            v = vg
            y = yg
            if v is None:
                v = x[mp[el[0,0]]]
            if y is None:
                y = x[mp[el[1,0]]]
            return y-np.log(norm.cdf(v))

        # grad f
        def grad_f(x):
            v = vg
            g = opt.spmatrix(0.0,[],[],(n,1))
            if v is None:
                v = x[mp[el[0,0]]]
            if type(el[0,0]) is cvxpy_scalar_var:
                g[mp[el[0,0]]] = -norm.pdf(v)/norm.cdf(v)
            if type(el[1,0]) is cvxpy_scalar_var:
                g[mp[el[1,0]]] = 1.
            return g

        # hess f
        def hess_f(x):
            v = vg
            h = opt.spmatrix(0.0,[],[],(n,n))
            if v is None:
                v = x[mp[el[0,0]]]
            if type(el[0,0]) is cvxpy_scalar_var:
                t = norm.pdf(v)
                s = norm.cdf(v)
                r = t/s
                h[mp[el[0,0]],mp[el[0,0]]] = r **2.+v*r
            return h

        # ind f
        def ind_f(x):
            return True

        # Return functions
        return f, grad_f, hess_f, ind_f

# Create instance
log_norm_cdf_hypo = cvxpy_log_norm_cdf_hypo()
