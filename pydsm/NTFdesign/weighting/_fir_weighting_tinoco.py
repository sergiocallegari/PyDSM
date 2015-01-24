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
import cvxpy_tinoco


def ntf_fir_from_digested(Qs, A, C, H_inf, **opts):
    """
    Synthesize FIR NTF from predigested specification

    Version for the cvxpy_tinoco modeler.
    """
    quiet = not opts['show_progress']
    order = np.size(Qs, 0)-1
    br = cvxpy_tinoco.variable(order, 1, name='br')
    b = cvxpy_tinoco.vstack((1, br))
    X = cvxpy_tinoco.variable(order, order, structure='symmetric', name='X')
    Qs = cvxpy_tinoco.matrix(Qs)
    target = cvxpy_tinoco.norm2(Qs*b)
    A = cvxpy_tinoco.matrix(A)
    B = cvxpy_tinoco.vstack((cvxpy_tinoco.zeros((order-1, 1)), 1.))
    C = cvxpy_tinoco.matrix(C)
    C = C+(cvxpy_tinoco.matrix(np.eye(order, order)[:, ::-1])*br).T
    D = cvxpy_tinoco.matrix(1.)
    M1 = A.T*X
    M2 = M1*B
    M = cvxpy_tinoco.vstack((
        cvxpy_tinoco.hstack((M1*A-X, M2, C.T)),
        cvxpy_tinoco.hstack((M2.T, B.T*X*B-H_inf**2, D)),
        cvxpy_tinoco.hstack((C, D, -1))
        ))
    constraint1 = cvxpy_tinoco.belongs(-M, cvxpy_tinoco.semidefinite_cone)
    constraint2 = cvxpy_tinoco.belongs(X, cvxpy_tinoco.semidefinite_cone)
    p = cvxpy_tinoco.program(cvxpy_tinoco.minimize(target),
                             [constraint1, constraint2])
    p.options.update(opts["tinoco_opts"])
    p.solve(quiet)
    return np.hstack((1, np.asarray(br.value.T)[0]))
