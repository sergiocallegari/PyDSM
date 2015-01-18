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
import cvxpy
from ...utilities import mdot
from ...delsig import padr


def ntf_fir_from_q0(q0, H_inf=1.5, **opts):
    """
    Synthesize FIR NTF from quadratic form expressing noise weighting.

    Version for the cvxpy modeler.
    """
    verbose = opts['show_progress']
    # Do the computation
    order = q0.shape[0]-1
    Q = la.toeplitz(q0)
    d, v = np.linalg.eigh(Q)
    if opts['fix_pos']:
        d = d/np.max(d)
        d[d < 0] = 0.
    qs = np.matrix(mdot(v, np.diag(np.sqrt(d)), np.linalg.inv(v)))
    br = cvxpy.Variable(order, 1, name='br')
    b = cvxpy.vstack(1, br)
    target = cvxpy.Minimize(cvxpy.norm2(qs*b))
    X = cvxpy.Semidef(order, name='X')
    A = np.matrix(np.eye(order, order, 1))
    B = np.vstack((np.zeros((order-1, 1)), 1.))
    C = br[::-1].T
    D = np.matrix(1.)
    M1 = A.T*X
    M2 = M1*B
    M = cvxpy.vstack(
        cvxpy.hstack(M1*A-X, M2, C.T),
        cvxpy.hstack(M2.T, B.T*X*B-H_inf**2, D),
        cvxpy.hstack(C, D, np.matrix(-1.))
        )
    Mconstr = -cvxpy.Semidef(order+2)
    constraints = [
        cvxpy.diag(Mconstr) == cvxpy.diag(M),
        cvxpy.upper_tri(Mconstr) == cvxpy.upper_tri(M),
    ]
    p = cvxpy.Problem(target, constraints)
    p.solve(solver=cvxpy.CVXOPT, verbose=verbose,
            kktsolver="chol", **opts['cvxopt_opts'])
    ntf_ir = np.hstack((1, np.asarray(br.value.T)[0]))
    return (np.roots(ntf_ir), np.zeros(order), 1.)


def ntf_hybrid_from_q0(q0, H_inf=1.5, poles=[], **opts):
    """
    Synthesize NTF from quadratic form expressing noise weighting and poles.

    Version for the cvxpy_tinoco modeler.
    """
    verbose = opts['show_progress']
    # Do the computation
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
    qs = np.matrix(mdot(v, np.diag(np.sqrt(d)), np.linalg.inv(v)))
    br = cvxpy.Variable(order, 1, name='br')
    b = cvxpy.vstack(1, br)
    target = cvxpy.Minimize(cvxpy.norm2(qs*b))
    X = cvxpy.Semidef(order, name='X')
    A = np.matrix(np.eye(order, order, 1))
    A[order-1] = -ar[::-1]
    B = np.vstack((np.zeros((order-1, 1)), 1.))
    C = br[::-1].T - ar[::-1].reshape((1, order))
    D = np.matrix(1.)
    M1 = A.T*X
    M2 = M1*B
    M = cvxpy.vstack(
        cvxpy.hstack(M1*A-X, M2, C.T),
        cvxpy.hstack(M2.T, B.T*X*B-H_inf**2, D),
        cvxpy.hstack(C, D, np.matrix(-1.))
        )
    Mconstr = -cvxpy.Semidef(order+2)
    constraints = [
        cvxpy.diag(Mconstr) == cvxpy.diag(M),
        cvxpy.upper_tri(Mconstr) == cvxpy.upper_tri(M),
    ]
    p = cvxpy.Problem(target, constraints)
    p.solve(solver=cvxpy.CVXOPT, verbose=verbose,
            kktsolver="chol", **opts['cvxopt_opts'])
    ntf_ir = np.hstack((1, np.asarray(br.value.T)[0]))
    return (np.roots(ntf_ir), poles, 1.)
