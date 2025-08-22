# -*- coding: utf-8 -*-

# Copyright (c) 2013–2024, Sergio Callegari
# All rights reserved.

# This file is part of PyDSM.

# PyDSM is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# PyDSM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with PyDSM.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import picos
import cvxopt


def ntf_fir_from_digested(Qs, A, C, H_inf, **opts):
    """
    Synthesize FIR NTF from predigested specification

    Version for the picos modeler.
    """
    verbose = 1 if opts.get('show_progress', True) else 0
    if 'maxiters' in opts['picos_opts']:
        opts['picos_opts']['max_iterations'] = opts['picos_opts']['maxiters']
        del opts['picos_opts']['maxiters']
    opts['picos_opts']['abs_ipm_opt_tol']=opts['picos_opts']['abstol']
    del opts['picos_opts']['abstol']
    opts['picos_opts']['rel_ipm_opt_tol']=opts['picos_opts']['reltol']
    del opts['picos_opts']['reltol']
    del opts['picos_opts']['feastol']
    # Do the computation
    order = int(np.size(Qs, 0)-1)
    Qs = cvxopt.matrix(Qs)
    p = picos.Problem(solver='cvxopt')
    power = picos.RealVariable('pow')
    br = picos.RealVariable('br', (order, 1))
    b = 1 // br
    X = picos.SymmetricVariable('X', order)
    A = cvxopt.matrix(A)
    B = cvxopt.matrix(np.vstack((np.zeros((order-1, 1)), 1.)))
    C = cvxopt.matrix(C)+br[::-1, :].T
    D = cvxopt.matrix(1.)
    M1 = A.T*X
    M2 = M1*B
    M = (((M1*A-X) & M2 & C.T) //
         (M2.T & (B.T*X*B-H_inf**2) & D) //
         (C & D & -1.))
    constraint1 = (M << 0)
    constraint2 = (X >> 0)
    p.set_objective('min', power)
    p.add_constraint(abs(Qs*b) <= power)
    p.add_constraint(constraint1)
    p.add_constraint(constraint2)
    p.solve(verbosity=verbose, **opts['picos_opts'])
    return np.hstack((1, np.asarray(br.value.T)[0]))
