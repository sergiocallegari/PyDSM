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
# along with PyDSM.  If not, see <https://www.gnu.org/licenses/>.

# This file includes code ported from the DELSIG Matlab toolbox
# (see https://www.mathworks.com/matlabcentral/fileexchange/19)
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

from __future__ import division

import numpy as np

from ._tf import evalTF


__all__ = ["rmsGain"]


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

    Examples
    --------

    >>> import numpy as np
    >>> import pydsm
    >>> v = np.sin(np.linspace(0, 20*np.pi, 1024))
    >>> print(pydsm.delsig.rms(v))
    0.7067614298686651

    Notes
    -----
    The discrete-time transer function H is evaluated in the frequency band
    (f1,f2).  Spanning of the bandwidth is linear. Frequencies are normalized
    in the [0, 0.5] interval.

    Warning: the result of the computation is normalized in the number of
    points used for the computation, but not on the frequency range.

    This computation could be practiced much more accurately and possibly
    faster, using algorithms for the numerical computation of integrals.
    """
    w = np.linspace(2*np.pi*f1, 2*np.pi*f2, N)
    return np.linalg.norm(evalTF(H, np.exp(1j*w))) / np.sqrt(N)
