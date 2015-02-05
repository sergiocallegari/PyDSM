# -*- coding: utf-8 -*-

# Copyright (c) 2015, Sergio Callegari
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

# This file includes portions of code ported from the "Frequency-Domain
# Min-Max Optimization for Delta-Sigma Modulators" Matlab toolbox
# (see http://www.mathworks.com/matlabcentral/fileexchange/36187)
# covered by the following copyright and permission notice
#
# Copyright (c) 2012 Masaaki Nagahara
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the distribution
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


from __future__ import division, print_function

import numpy as np
import cvxpy_tinoco


def ntf_fir_from_digested(order, osrs, H_inf, f0s, zf, **opts):
    """
    Synthesize FIR NTF with minmax approach from predigested specification

    Version for the cvxpy_tinoco modeler.
    """
    quiet = not opts['show_progress']

    # State space representation of NTF
    A = cvxpy_tinoco.matrix(np.eye(order, order, 1))
    B = cvxpy_tinoco.matrix(np.vstack((np.zeros((order-1, 1)), 1.)))
    # C contains the NTF coefficients
    D = cvxpy_tinoco.matrix([[1]])

    # Set up the problem
    bands = len(f0s)
    c = cvxpy_tinoco.variable(1, order)
    F = []
    gg = cvxpy_tinoco.arrays.cvxpy_var(bands, 1)

    for idx in range(bands):
        f0 = f0s[idx]
        osr = osrs[idx]
        omega0 = 2*f0*np.pi
        Omega = 1./osr*np.pi
        P = cvxpy_tinoco.variable(order, order, 'symmetric')
        Q = cvxpy_tinoco.variable(order, order, 'symmetric')
        if f0 == 0:
            # Lowpass modulator
            M1 = A.T*P*A+Q*A+A.T*Q-P-2*Q*np.cos(Omega)
            M2 = A.T*P*B + Q*B
            M3 = B.T*P*B - gg[idx, 0]
            M = cvxpy_tinoco.vstack(
                (cvxpy_tinoco.hstack((M1, M2, c.T)),
                 cvxpy_tinoco.hstack((M2.T, M3, D)),
                 cvxpy_tinoco.hstack((c, D, -1))))
            F += [cvxpy_tinoco.belongs(Q, cvxpy_tinoco.semidefinite_cone)]
            F += [cvxpy_tinoco.belongs(-M, cvxpy_tinoco.semidefinite_cone)]
            if zf:
                # Force a zero at DC
                F += [cvxpy_tinoco.equals(cvxpy_tinoco.sum(c), -1)]
        else:
            # Bandpass modulator
            M1r = (A.T*P*A + Q*A*np.cos(omega0) + A.T*Q*np.cos(omega0) -
                   P - 2*Q*np.cos(Omega))
            M2r = A.T*P*B + Q*B*np.cos(omega0)
            M3r = B.T*P*B - gg[idx, 0]
            M1i = A.T*Q*np.sin(omega0) - Q*A*np.sin(omega0)
            M21i = -Q*B*np.sin(omega0)
            M22i = B.T*Q*np.sin(omega0)
            Mr = cvxpy_tinoco.vstack(
                (cvxpy_tinoco.hstack((M1r, M2r, c.T)),
                 cvxpy_tinoco.hstack((M2r.T, M3r, D)),
                 cvxpy_tinoco.hstack((c, D, -1))))
            Mi = cvxpy_tinoco.vstack(
                (cvxpy_tinoco.hstack((M1i, M21i,
                                      cvxpy_tinoco.zeros((order, 1)))),
                 cvxpy_tinoco.hstack((M22i, 0, 0)),
                 cvxpy_tinoco.hstack((cvxpy_tinoco.zeros((1, order)), 0, 0))))
            M = cvxpy_tinoco.vstack(
                (cvxpy_tinoco.hstack((Mr, Mi)),
                 cvxpy_tinoco.hstack((-Mi, Mr))))
            F += [cvxpy_tinoco.belongs(Q, cvxpy_tinoco.semidefinite_cone)]
            F += [cvxpy_tinoco.belongs(-M, cvxpy_tinoco.semidefinite_cone)]
            if zf:
                # Force a zero at z=np.exp(1j*omega0)
                nn = np.arange(order).reshape((order, 1))
                vr = cvxpy_tinoco.matrix(np.cos(omega0*nn))
                vi = cvxpy_tinoco.matrix(np.sin(omega0*nn))
                vn = cvxpy_tinoco.matrix(
                    [-np.cos(omega0*order), -np.sin(omega0*order)])
                F += [cvxpy_tinoco.equals(c*cvxpy_tinoco.hstack((vr, vi)), vn)]
    if H_inf < np.inf:
        # Enforce the Lee constraint
        R = cvxpy_tinoco.variable(order, order, 'symmetric')
        F += [cvxpy_tinoco.belongs(R, cvxpy_tinoco.semidefinite_cone)]
        MM = cvxpy_tinoco.vstack(
            (cvxpy_tinoco.hstack((A.T*R*A-R, A.T*R*B, c.T)),
             cvxpy_tinoco.hstack((B.T*R*A, -H_inf**2+B.T*R*B, D)),
             cvxpy_tinoco.hstack((c, D, -1))))
        F += [cvxpy_tinoco.belongs(-MM, cvxpy_tinoco.semidefinite_cone)]
    p = cvxpy_tinoco.program(cvxpy_tinoco.minimize(cvxpy_tinoco.max(gg)), F)
    p.options.update(opts["tinoco_opts"])
    p.solve(quiet)
    return np.hstack((1, np.asarray(c.value)[0, ::-1]))
