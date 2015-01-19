# -*- coding: utf-8 -*-

# Copyright (c) 2013, Sergio Callegari
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
# along with PyDSM.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import scipy.linalg as la
import picos
import cvxopt
from ...utilities import mdot
from ...delsig import padr


def ntf_fir_from_q0(q0, H_inf=1.5, **opts):
    """
    Synthesize FIR NTF from quadratic form expressing noise weighting.

    Version for the cvxpy_tinoco modeler.
    """
    order = q0.shape[0]-1
    Q = la.toeplitz(q0)
    d, v = np.linalg.eigh(Q)
    if opts['fix_pos']:
        d = d/np.max(d)
        d[d < 0] = 0.
    Qs = mdot(v, np.diag(np.sqrt(d)), np.linalg.inv(v))
    A = np.eye(order, order, 1)
    C = np.zeros((1, order))
    ntf_ir = ntf_fir_from_digested(Qs, A, C, H_inf=1.5, **opts)
    return (np.roots(ntf_ir), np.zeros(order), 1.)


def ntf_hybrid_from_q0(q0, H_inf=1.5, poles=[], **opts):
    """
    Synthesize NTF from quadratic form expressing noise weighting and poles.

    Version for the cvxpy_tinoco modeler.
    """
    order = q0.shape[0]-1
    poles = np.asarray(poles).reshape(-1)
    if poles.shape[0] > order:
        raise ValueError('Too many poles provided')
    poles = padr(poles, order, 0)
    # Get denominator coefficients from a_1 to a_order (a_0 is 1)
    ar = np.poly(poles)[1:].real
    Q = la.toeplitz(q0)
    d, v = np.linalg.eigh(Q)
    if opts['fix_pos']:
        d = d/np.max(d)
        d[d < 0] = 0.
    Qs = mdot(v, np.diag(np.sqrt(d)), np.linalg.inv(v))
    A = np.eye(order, order, 1)
    A[order-1] = -ar[::-1]
    C = -ar[::-1].reshape((1, order))
    ntf_ir = ntf_fir_from_digested(Qs, A, C, H_inf=1.5, **opts)
    return (np.roots(ntf_ir), poles, 1.)


def ntf_fir_from_digested(Qs, A, C, H_inf=1.5, **opts):
    """
    Synthesize FIR NTF from predigested specification

    Version for the cvxpy modeler.
    """
    verbose = 1 if opts.get('show_progress', True) else 0
    cvxopt_opts = opts['cvxopt_opts']
    if 'maxiters' in cvxopt_opts:
        cvxopt_opts['maxit'] = cvxopt_opts['maxiters']
        del cvxopt_opts['maxiters']
    # Do the computation
    order = np.size(Qs, 0)-1
    Qs = cvxopt.matrix(Qs)
    p = picos.Problem(solver='cvxopt')
    power = p.add_variable('pow')
    br = p.add_variable('br', (order, 1))
    b = 1 // br
    X = p.add_variable('X', (order, order), vtype='symmetric')
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
    p.add_constraint(abs(Qs*b) < power)
    p.add_constraint(constraint1)
    p.add_constraint(constraint2)
    p.set_options(**opts['cvxopt_opts'])
    p.solve(verbose=verbose)
    return np.hstack((1, np.asarray(br.value.T)[0]))
