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
import cvxpy_tinoco
from ...utilities import mdot
from ...delsig import padr


def ntf_fir_from_q0(q0, H_inf=1.5, **opts):
    """
    Synthesize FIR NTF from quadratic form expressing noise weighting.

    Version for the cvxpy_tinoco modeler.
    """
    quiet = not opts['show_progress']
    # Do the computation
    order = q0.shape[0]-1
    Q = la.toeplitz(q0)
    d, v = np.linalg.eigh(Q)
    if opts['fix_pos']:
        d = d/np.max(d)
        d[d < 0] = 0.
    qs = cvxpy_tinoco.matrix(mdot(v, np.diag(np.sqrt(d)), np.linalg.inv(v)))
    br = cvxpy_tinoco.variable(order, 1, name='br')
    b = cvxpy_tinoco.vstack((1, br))
    target = cvxpy_tinoco.norm2(qs*b)
    X = cvxpy_tinoco.variable(order, order, structure='symmetric', name='X')
    A = cvxpy_tinoco.matrix(np.eye(order, order, 1))
    B = cvxpy_tinoco.vstack((cvxpy_tinoco.zeros((order-1, 1)), 1.))
    C = (cvxpy_tinoco.matrix(np.eye(order, order)[:, ::-1])*br).T
    D = cvxpy_tinoco.matrix(1.)
    M1 = A.T*X
    M2 = M1*B
    M = cvxpy_tinoco.vstack((
        cvxpy_tinoco.hstack((M1*A-X, M2, C.T)),
        cvxpy_tinoco.hstack((M2.T, B.T*X*B-H_inf**2, D)),
        cvxpy_tinoco.hstack((C, D, cvxpy_tinoco.matrix(-1.)))
        ))
    constraint1 = cvxpy_tinoco.belongs(-M, cvxpy_tinoco.semidefinite_cone)
    constraint2 = cvxpy_tinoco.belongs(X, cvxpy_tinoco.semidefinite_cone)
    p = cvxpy_tinoco.program(cvxpy_tinoco.minimize(target),
                             [constraint1, constraint2])
    p.options.update(opts["cvxopt_opts"])
    p.solve(quiet)
    ntf_ir = np.hstack((1, np.asarray(br.value.T)[0]))
    return (np.roots(ntf_ir), np.zeros(order), 1.)


def ntf_hybrid_from_q0(q0, H_inf=1.5, poles=[], **opts):
    """
    Synthesize NTF from quadratic form expressing noise weighting and poles.

    Version for the cvxpy_tinoco modeler.
    """
    quiet = not opts['show_progress']
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
    qs = cvxpy_tinoco.matrix(mdot(v, np.diag(np.sqrt(d)), np.linalg.inv(v)))
    br = cvxpy_tinoco.variable(order, 1, name='br')
    b = cvxpy_tinoco.vstack((1, br))
    target = cvxpy_tinoco.norm2(qs*b)
    X = cvxpy_tinoco.variable(order, order, structure='symmetric', name='X')
    A = cvxpy_tinoco.matrix(np.eye(order, order, 1))
    A[order-1] = -ar[::-1]
    B = cvxpy_tinoco.vstack((cvxpy_tinoco.zeros((order-1, 1)), 1.))
    C = (cvxpy_tinoco.matrix(np.eye(order, order)[:, ::-1])*br).T
    C = C-cvxpy_tinoco.matrix(ar[::-1])
    D = cvxpy_tinoco.matrix(1.)
    M1 = A.T*X
    M2 = M1*B
    M = cvxpy_tinoco.vstack((
        cvxpy_tinoco.hstack((M1*A-X, M2, C.T)),
        cvxpy_tinoco.hstack((M2.T, B.T*X*B-H_inf**2, D)),
        cvxpy_tinoco.hstack((C, D, cvxpy_tinoco.matrix(-1.)))
        ))
    constraint1 = cvxpy_tinoco.belongs(-M, cvxpy_tinoco.semidefinite_cone)
    constraint2 = cvxpy_tinoco.belongs(X, cvxpy_tinoco.semidefinite_cone)
    p = cvxpy_tinoco.program(cvxpy_tinoco.minimize(target),
                             [constraint1, constraint2])
    p.options.update(opts["cvxopt_opts"])
    p.solve(quiet)
    ntf_ir = np.hstack((1, np.asarray(br.value.T)[0]))
    return (np.roots(ntf_ir), poles, 1.)
