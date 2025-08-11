# -*- coding: utf-8 -*-

# Copyright (c) 2015–2024, Sergio Callegari
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
import cvxpy


def ntf_fir_from_digested(order, osrs, H_inf, f0s, zf, **opts):
    """
    Synthesize FIR NTF with minmax approach from predigested specification

    Version for the cvxpy modeler.
    """
    verbose = opts['show_progress']
    if opts['cvxpy_opts']['solver'] == 'cvxopt':
        opts['cvxpy_opts']['solver'] = cvxpy.CVXOPT
    elif opts['cvxpy_opts']['solver'] == 'scs':
        opts['cvxpy_opts']['solver'] = cvxpy.SCS

    # State space representation of NTF
    A = np.eye(order, order, 1)
    B = np.vstack((np.zeros((order-1, 1)), 1.))
    # C contains the NTF coefficients
    D = np.array([[1.]])

    # Set up the problem
    bands = len(f0s)
    c = cvxpy.Variable((1, order))
    F = []
    gg = cvxpy.Variable((bands, 1))

    for idx in range(bands):
        f0 = f0s[idx]
        osr = osrs[idx]
        omega0 = 2*f0*np.pi
        Omega = 1./osr*np.pi
        P = cvxpy.Variable((order, order), symmetric=True)
        Q = cvxpy.Variable((order, order), PSD=True)
        if f0 == 0:
            # Lowpass modulator
            M1 = (A.T @ P @ A + Q @ A +
                  A.T @ Q - P - 2*Q*np.cos(Omega))
            M2 = A.T @ P @ B + Q @ B
            M3 = B.T @ P @ B - gg[idx, 0]
            M = cvxpy.bmat([[M1, M2, c.T],
                            [M2.T, M3, D],
                            [c, D, np.array([[-1]])]])
            F += [M << 0]
            if zf:
                # Force a zero at DC
                F += [cvxpy.sum_entries(c) == -1]
        else:
            # Bandpass modulator
            M1r = (A.T @ P @ A + Q @ A*np.cos(omega0) +
                   A.T @ Q*np.cos(omega0) -
                   P - 2*Q*np.cos(Omega))
            M2r = A.T @ P @ B + Q @ B*np.cos(omega0)
            M3r = B.T @ P @ B - gg[idx, 0]
            M1i = A.T @ Q*np.sin(omega0) - Q @ A*np.sin(omega0)
            M21i = -Q @ B*np.sin(omega0)
            M22i = B.T @ Q*np.sin(omega0)
            Mr = cvxpy.bmat([[M1r, M2r, c.T],
                             [M2r.T, M3r, D],
                             [c, D, np.array([[-1]])]])
            Mi = cvxpy.bmat([[M1i, M21i, np.zeros((order, 1))],
                             [M22i, np.array([[0]]), np.array([[0]])],
                             [np.zeros((1, order)), np.array([[0]]),
                              np.array([[0]])]])
            M = cvxpy.bmat([[Mr, Mi],
                            [-Mi, Mr]])
            F += [M << 0]
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
        R = cvxpy.Variable((order, order), PSD=True)
        MM = cvxpy.bmat([[A.T @ R @ A - R, A.T @ R @ B, c.T],
                         [B.T @ R @ A, -H_inf**2 + B.T @ R @ B, D],
                         [c, D, np.array([[-1]])]])
        F += [MM << 0]
    target = cvxpy.Minimize(cvxpy.max(gg))
    p = cvxpy.Problem(target, F)
    p.solve(verbose=verbose, **opts['cvxpy_opts'])
    return np.hstack((1, np.asarray(c.value)[0, ::-1]))
