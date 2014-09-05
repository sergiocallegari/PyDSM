# -*- coding: utf-8 -*-

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
Transfer function management a la DELSIG.
=========================================
"""

import numpy as np

__all__=["evalTF", "evalRPoly", "rmsGain"]

def evalTF(tf, x):
    """
    Evaluates a transfer function.

    This function can be used either for discrete time or continuous time
    transfer functions (see the notes).

    Parameters
    ----------
    tf : tuple
        transfer function in zpk or ba form
    x : complex or array_like of complex
        value or vector of values where the tf is to be evaluated

    Returns
    -------
    y : ndarray
        value(s) of transfer function at the given complex values.

    Notes
    -----
    Parameter x corresponds to 's' or 'z' in CT or DT transfer functions
    respectively. Thus it should be 1j*omega or exp(1j*omega*T).

    With respect to the tf parameter, zpk form is a triple containing a list
    of zeros, a list of poles and a scalar gain. The ba form (also called tf
    form in the scipy documentation) is a couple containing a list of the
    numerator coefficients and a list of the denominator coefficients. The
    coefficients are sorted from the higher power of 's' or 'z' to the lower,
    so that the last coefficient is in fact the constant term in the
    numerator/denominator polynomial.
    """
    if len(tf) == 3:
        retval = np.asarray(evalRPoly(tf[0], x, tf[2]))/evalRPoly(tf[1], x, 1)
    elif len(tf) == 2:
        retval = np.asarray(np.polyval(tf[0], x))/np.polyval(tf[1],x)
    if retval.size == 1:
        return retval.item()
    return retval

def evalRPoly(roots, x, k=1):
    """
    Compute the value of a polynomial that is given in terms of its roots.
    Roots at infinity are removed before the computation.

    Parameters
    ----------
    roots : array_like
        roots of polynomial
    x : complex or array_like of complex
        complex value or vector of complex values where tf is to be
        evaluated
    k -> real, optional
        gain, defaults to 1.

    Returns
    -------
    y : ndarray
        value of polynomial at the given complex values.
    """
    y = k*np.ones_like(x)
    roots = np.asarray(roots)
    # Remove roots at infinity
    roots = roots[np.logical_not(np.isinf(roots))]
    for i in xrange(roots.size):
        y = y*(x-roots[i])
    return y

def rmsGain(H, f1, f2, N=100):
    """
    Compute the root-mean-square gain of a DT transfer function.

    The computation is practiced on an assigned frequency interval.
    By root-mean-square, it is meant that the magnitude response of the
    transfer function is squared and averaged over the band of interest.
    Then, the square root of this averaged value is returned.

    Parameters
    ----------
    H : tuple
        transfer function either in (z,p,k) or (n,d) form
    f1 : real
        lower bound of frequency band on which the transfer function
        is evaluated
    f2 : real
        upper bound of frequency band on which the transfer function
        is evaluated
    N : int, optional
        number of points where the transfer function is evaluated in the
        interval. Defaults to 100.

    Returns
    -------
    rms : real
        rms value of the discrete time transfer function.

    Notes
    -----
    The discrete-time transer function H is evaluated in the frequency band
    (f1,f2).  Spanning of the bandwidth is linear. Frequencies are normalized
    in the [0,0.5] interval.

    Warning: the result of the computation is normalized in the number of
    points used for the computation, but not on the frequency range.

    This computation could be practiced much more accurately and possibly
    faster, using algorithms for the numerical computation of integrals.
    """
    w = np.linspace(2*np.pi*f1, 2*np.pi*f2, N)
    return np.linalg.norm(evalTF(H, np.exp(1j*w))) / np.sqrt(N)
