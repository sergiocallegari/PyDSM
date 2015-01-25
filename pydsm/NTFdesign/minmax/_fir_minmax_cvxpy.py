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

from __future__ import division, print_function

import numpy as np
import cvxpy


def ntf_fir_from_digested(order, osrs, H_inf, f0s, zf, **opts):
    """
    Synthesize FIR NTF with minmax approach from predigested specification

    Version for the cvxpy_tinoco modeler.
    """
    verbose = opts['show_progress']
    if opts['cvxpy_opts']['solver'] == 'cvxopt':
        opts['cvxpy_opts']['solver'] = cvxpy.CVXOPT
    elif opts['cvxpy_opts']['solver'] == 'scs':
        opts['cvxpy_opts']['solver'] = cvxpy.SCS

    # State space representation of NTF
    A = np.matrix(np.eye(order, order, 1))
    B = np.matrix(np.vstack((np.zeros((order-1, 1)), 1.)))
    # C contains the NTF coefficients
    D = np.matrix(1)

    # Set up the problem
    bands = len(f0s)
    c = cvxpy.Variable(1, order)
    F = []
    gg = cvxpy.Variable(bands, 1)

    for idx in range(bands):
        f0 = f0s[idx]
        osr = osrs[idx]
        omega0 = 2*f0*np.pi
        Omega = 1./osr*np.pi
        P = cvxpy.Variable(order, order)
        F += [cvxpy.upper_tri(P) == cvxpy.upper_tri(P.T)]
        Q = cvxpy.Semidef(order)
        if f0 == 0:
            # Lowpass modulator
            M1 = A.T*P*A+Q*A+A.T*Q-P-2*Q*np.cos(Omega)
            M2 = A.T*P*B + Q*B
            M3 = B.T*P*B - gg[idx, 0]
            M = cvxpy.vstack(
                cvxpy.hstack(M1, M2, c.T),
                cvxpy.hstack(M2.T, M3, D),
                cvxpy.hstack(c, D, -1))
            Mconstr = -cvxpy.Semidef(order+2)
            F += [cvxpy.diag(Mconstr) == cvxpy.diag(M),
                  cvxpy.upper_tri(Mconstr) == cvxpy.upper_tri(M)]
            if zf:
                # Force a zero at DC
                F += [cvxpy.sum_entries(c) == -1]
        else:
            # Bandpass modulator
            M1r = (A.T*P*A + Q*A*np.cos(omega0) + A.T*Q*np.cos(omega0) -
                   P - 2*Q*np.cos(Omega))
            M2r = A.T*P*B + Q*B*np.cos(omega0)
            M3r = B.T*P*B - gg[idx, 0]
            M1i = A.T*Q*np.sin(omega0) - Q*A*np.sin(omega0)
            M21i = -Q*B*np.sin(omega0)
            M22i = B.T*Q*np.sin(omega0)
            Mr = cvxpy.vstack(
                cvxpy.hstack(M1r, M2r, c.T),
                cvxpy.hstack(M2r.T, M3r, D),
                cvxpy.hstack(c, D, -1))
            Mi = cvxpy.vstack(
                cvxpy.hstack(M1i, M21i, np.zeros((order, 1))),
                cvxpy.hstack(M22i, 0, 0),
                cvxpy.hstack(np.zeros((1, order)), 0, 0))
            M = cvxpy.vstack(
                cvxpy.hstack(Mr, Mi),
                cvxpy.hstack(-Mi, Mr))
            Mconstr = -cvxpy.Semidef(2*(order+2))
            F += [cvxpy.diag(Mconstr) == cvxpy.diag(M),
                  cvxpy.upper_tri(Mconstr) == cvxpy.upper_tri(M)]
            if zf:
                # Force a zero at z=np.exp(1j*omega0)
                nn = np.arange(order).reshape((order, 1))
                vr = np.matrix(np.cos(omega0*nn))
                vi = np.matrix(np.sin(omega0*nn))
                vn = np.matrix(
                    [-np.cos(omega0*order), -np.sin(omega0*order)])
                F += [c*cvxpy.hstack(vr, vi) == vn]
    if H_inf < np.inf:
        # Enforce the Lee constraint
        R = cvxpy.Semidef(order)
        MM = cvxpy.vstack(
            cvxpy.hstack(A.T*R*A-R, A.T*R*B, c.T),
            cvxpy.hstack(B.T*R*A, -H_inf**2+B.T*R*B, D),
            cvxpy.hstack(c, D, -1))
        Mconstr = -cvxpy.Semidef(order+2)
        F += [cvxpy.diag(Mconstr) == cvxpy.diag(MM),
              cvxpy.upper_tri(Mconstr) == cvxpy.upper_tri(MM)]
    target = cvxpy.Minimize(cvxpy.max_entries(gg))
    p = cvxpy.Problem(target, F)
    p.solve(verbose=verbose, **opts['cvxpy_opts'])
    return np.hstack((1, np.asarray(c.value)[0, ::-1]))
