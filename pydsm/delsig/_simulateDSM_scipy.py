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

"""
Slow simulator for a generic delta sigma modulator
==================================================
"""

import numpy as np
import scipy as sp
__import__('scipy.signal')
from warnings import warn
from ..errors import PyDsmWarning, PyDsmError

__all__=["ds_quantize"]

def simulateDSM(u, arg2, nlev=2, x0=0,
                store_xn=False, store_xmax=False, store_y=False):

    warn('Running the slow version of simulateDSM.', PyDsmWarning)

    # Make sure that nlev is an array
    nlev=np.asarray(nlev).reshape(1)

    # Make sure that input is a matrix
    u = np.asarray(u)
    if u.ndim == 1:
        u = u.reshape(1, -1)

    nu = u.shape[0]
    nq = np.size(nlev)

    if type(arg2)==tuple and len(arg2)==3:
        # Assume ntf in zpk form
        (ntf_z, ntf_p, ntf_k) = arg2
        form = 2
        order = len(ntf_z)
    else:
        # Assume ABCD form
        ABCD = np.asarray(arg2)
        if ABCD.shape[1] == nu+ABCD.shape[0]:
            # ABCD dimensions OK
            form = 1
            order = ABCD.shape[0]-nq
        else:
            raise PyDsmError('Incorrect modulator specification')

    # Assure that the state is a column vector
    if np.isscalar(x0) and x0 == 0:
        x0 = np.zeros((order, 1))
    else:
        x0 = np.array(x0, dtype=float).reshape(-1, 1)

    if form == 1:
        A = ABCD[0:order, 0:order]
        B = ABCD[0:order, order:order+nu+nq]
        C = ABCD[order:order+nq, 0:order]
        D1 = ABCD[order:order+nq, order:order+nu]
    else:
        # Seek a realization of -1/H
        A, B2, C, D2 = sp.signal.zpk2ss(ntf_p, ntf_z, -1)
        # Transform the realization so that C = [1 0 0 ...]
        Sinv = sp.linalg.orth(np.hstack((np.transpose(C), np.eye(order))))/ \
            np.linalg.norm(C)
        S = sp.linalg.inv(Sinv)
        C = np.dot(C, Sinv)
        if C[0, 0] < 0:
            S = -S
            Sinv = -Sinv
        A = np.dot(np.dot(S, A), Sinv)
        B2 = np.dot(S, B2)
        C = np.hstack(([[1]], np.zeros((1,order-1))))
        # C=C*Sinv;
        # D2 = 0;
        # !!!! Assume stf=1
        B1 = -B2
        D1 = 1
        B = np.hstack((B1, B2))

    N = u.shape[1]
    v = np.empty((nq, N))
    if store_y:
        # Need to store the quantizer input
        y = np.empty((nq, N))
    else:
        y = np.empty((0,0))
    if store_xn:
        # Need to store the state information
        xn = np.empty((order, N))
    if store_xmax:
        # Need to keep track of the state maxima
        xmax = np.abs(x0)
    else:
        xmax = np.empty(0)

    for i in xrange(N):
        # I guess the coefficients in A, B, C, D should be real...
        y0 = np.real(np.dot(C, x0) + np.dot(D1, u[:, i]))
        if store_y:
            y[:, i] = y0
        v[:, i] = ds_quantize(y0, nlev)
        x0 = np.dot(A, x0) + np.dot(B, np.vstack((u[:, i], v[:, i])))
        if store_xn:
            # Save the next state
            xn[:, i] = x0
        if store_xmax:
            # Keep track of the state maxima
            xmax = np.max((np.abs(x0), xmax), 0)
    if not store_xn:
        xn = x0
    return v.squeeze(), xn.squeeze(), xmax, y.squeeze()

def ds_quantize(y, n):
    """Quantize a signal according to a given number of levels.

    Parameters
    ----------
    y : real or array of reals
        signal to be quantized (1 sample!). A column vector with more than
        1 row if there are multiple quantizers.
    n : int or vector of ints
        number of quantization levels. Can be a vector to specify multiple
        quantizers, in this case, y must have as many rows as the entries in
        n

    Returns
    -------
    z : real or ndarray
        quantized signal (1 sample!). A column vector with
        more than 1 row if there are multiple quantizers.

    Notes
    -----
    y is quantized to:

    * an odd integer in [-n+1, n-1], if n is even, or
    * an even integer in [-n, n], if n is odd.

    This definition gives the same step height for both mid-riser and
    mid-tread quantizers.
    """
    v=np.empty_like(y)
    for qi in xrange(np.size(n)):
        if np.remainder(n[qi], 2) == 0:
            v[qi] = 2*np.floor(0.5*y[qi])+1
        else:
            v[qi] = 2*np.floor(0.5*(y[qi]+1))
        L = n[qi]-1
        v[qi, 0]=np.max((np.min((v[qi, 0], L)), -L))
    return v
