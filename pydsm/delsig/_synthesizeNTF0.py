# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
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

# This file includes code ported from the DELSIG Matlab toolbox
# (see http://www.mathworks.com/matlabcentral/fileexchange/19)
# covered by the following copyright and permission notice
#
# Copyright (c) 2009 Richard Schreier
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the distribution
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import numpy as np
from warnings import warn
from ..exceptions import PyDsmApproximationWarning
from ._tf import evalTF
from ..relab import cplxpair
from ._ds import ds_optzeros


def synthesizeNTF0(order, osr, opt, H_inf, f0):
    # Determine the zeros.
    if f0 != 0:
        # Bandpass design-- halve the order temporarily.
        order = order/2
        dw = np.pi/(2*osr)
    else:
        dw = np.pi/osr

    if opt.ndim == 0:
        # opt is a number
        if opt == 0:
            z = np.zeros(order)
        else:
            z = dw*ds_optzeros(order, opt)
        if z.size == 0:
            raise RuntimeError('Cannot synthesize NTF zeros')
        if f0 != 0:
            # Bandpass design-- shift and replicate the zeros.
            order = order*2
            z = z + 2*np.pi*f0
            z = np.vstack((z, -z)).transpose().flatten()
        z = np.exp(1j*z)
    else:
        z = opt

    p = np.zeros(order)
    k = 1
    itn_limit = 100
    fprev = 0

    if f0 == 0:
        # Lowpass design
        HinfLimit = 2**order
        # !!! The limit is actually lower for opt=1 and low OSR
        if H_inf >= HinfLimit:
            warn('Unable to achieve specified H_inf, '
                 'setting all NTF poles to zero',
                 PyDsmApproximationWarning)
            p = np.zeros(order)
        else:
            x = 0.3**(order-1)   # starting guess
            for itn in xrange(1, itn_limit+1):
                me2 = -0.5*(x**(2./order))
                w = (2*np.arange(1, order+1)+1)*np.pi/order
                mb2 = 1+me2*np.exp(1j*w)
                p = mb2 - np.sqrt(mb2**2-1)
                # Reflect poles to be inside the unit circle
                out = abs(p) > 1
                p[out] = 1/p[out]
                # The following is not exactly what delsig does.
                # We do not have an identical cplxpair
                p = cplxpair(p)
                f = np.real(evalTF((z, p, k), -1))-H_inf
                if itn == 1:
                    delta_x = -f/100
                else:
                    delta_x = -f*delta_x/(f-fprev)
                xplus = x+delta_x
                if xplus > 0:
                    x = xplus
                else:
                    x = x*0.1
                fprev = f
                if abs(f) < 1e-10 or abs(delta_x) < 1e-10:
                    break
                if x > 1e6:
                    warn('Unable to achieve specified Hinf, '
                         'setting all NTF poles to zero',
                         PyDsmApproximationWarning)
                    p = np.zeros(order)
                    break
                if itn == itn_limit:
                    warn('Iteration limit exceeded',
                         PyDsmApproximationWarning)
    else:
        # Bandpass design
        x = 0.3**(order/2-1)   # starting guess (not very good for f0~0)
        if f0 > 0.25:
            z_inf = 1.
        else:
            z_inf = -1.
        c2pif0 = np.cos(2*np.pi*f0)
        for itn in xrange(1, itn_limit+1):
            e2 = 0.5*x**(2./order)
            w = (2*np.arange(order)+1)*np.pi/order
            mb2 = c2pif0 + e2*np.exp(1j*w)
            p = mb2 - np.sqrt(mb2**2-1)
            # Reflect poles to be inside the unit circle
            out = abs(p) > 1
            p[out] = 1/p[out]
            # The following is not exactly what delsig does.
            p = cplxpair(p)
            f = np.real(evalTF((z, p, k), z_inf))-H_inf
            if itn == 1:
                delta_x = -f/100
            else:
                delta_x = -f*delta_x/(f-fprev)
            xplus = x+delta_x
            if xplus > 0:
                x = xplus
            else:
                x = x*0.1
            fprev = f
            if abs(f) < 1e-10 or abs(delta_x) < 1e-10:
                break
            if x > 1e6:
                warn('Unable to achieve specified Hinf, '
                     'Setting all NTF poles to zero',
                     PyDsmApproximationWarning)
                p = np.zeros(order)
                break
            if itn == itn_limit:
                warn('Iteration limit exceeded', PyDsmApproximationWarning)

    z = cplxpair(z)
    return z, p, k
