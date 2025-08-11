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
import picos
import cvxopt


def ntf_fir_from_digested(order, osrs, H_inf, f0s, zf, **opts):
    """
    Synthesize FIR NTF with minmax approach from predigested specification

    Version for the PICOS modeler.
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

    # State space representation of NTF
    A = cvxopt.matrix(np.eye(order, order, 1))
    B = cvxopt.matrix(np.vstack((np.zeros((order-1, 1)), 1.)))
    # C contains the NTF coefficients
    D = cvxopt.matrix(1.)

    p = picos.Problem()
    # Set up the problem
    bands = len(f0s)
    c = picos.RealVariable('c', (1, order))
    gg = picos.RealVariable('gg', bands)
    #gg = p.add_variable('gg', bands)

    for idx in range(bands):
        f0 = f0s[idx]
        osr = osrs[idx]
        omega0 = 2*f0*np.pi
        Omega = 1./osr*np.pi
        P = picos.SymmetricVariable(f'P[{idx}]', order)
        Q = picos.SymmetricVariable(f'Q[{idx}]', order)
        if f0 == 0:
            # Lowpass modulator
            M1 = A.T*P*A+Q*A+A.T*Q-P-2*Q*np.cos(Omega)
            M2 = A.T*P*B + Q*B
            M3 = B.T*P*B - gg[idx, 0]
            M = ((M1 & M2 & c.T) //
                 (M2.T & M3 & D) //
                 (c & D & -1.0))
            p.add_constraint(Q >> 0)
            p.add_constraint(M << 0)
            if zf:
                # Force a zero at DC
                p.add_constraint(picos.sum(c[i] for i in range(order)) == -1)
        else:
            # Bandpass modulator
            M1r = (A.T*P*A + Q*A*np.cos(omega0) + A.T*Q*np.cos(omega0) -
                   P - 2*Q*np.cos(Omega))
            M2r = A.T*P*B + Q*B*np.cos(omega0)
            M3r = B.T*P*B - gg[idx]
            M1i = A.T*Q*np.sin(omega0) - Q*A*np.sin(omega0)
            M21i = -Q*B*np.sin(omega0)
            M22i = B.T*Q*np.sin(omega0)
            Mr = ((M1r & M2r & c.T) //
                  (M2r.T & M3r & D) //
                  (c & D & -1.0))
            Mi = ((M1i & M21i & cvxopt.matrix(np.zeros((order, 1)))) //
                  (M22i & 0.0 & 0.0) //
                  (cvxopt.matrix(np.zeros((1, order+2)))))
            M = ((Mr & Mi) //
                 (-Mi & Mr))
            p.add_constraint(Q >> 0.0)
            p.add_constraint(M << 0.0)
            if zf:
                # Force a zero at z=np.exp(1j*omega0)
                nn = np.arange(order).reshape((order, 1))
                vr = cvxopt.matrix(np.cos(omega0*nn))
                vi = cvxopt.matrix(np.sin(omega0*nn))
                vn = cvxopt.matrix(
                    [-np.cos(omega0*order), -np.sin(omega0*order)])
                p.add_constraint(c*(vr & vi) == vn)
    if H_inf < np.inf:
        # Enforce the Lee constraint
        R = picos.SymmetricVariable('R', order)
        p.add_constraint(R >> 0)
        MM = ((A.T*R*A-R & A.T*R*B & c.T) //
              (B.T*R*A & -H_inf**2+B.T*R*B & D) //
              (c & D & -1))
        p.add_constraint(MM << 0)
    mx = picos.RealVariable('mx')
    p.add_constraint(picos.sum_k_largest(gg, 1) <= mx)
    p.set_objective('min', mx)
    p.solve(verbosity=verbose, **opts['picos_opts'])
    return np.hstack((1, np.asarray(c.value)[0, ::-1]))
