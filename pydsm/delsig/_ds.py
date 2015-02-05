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
DELSIG helpers routines
=======================
"""

import numpy as np
from math import ceil, sqrt
from ._padding import padl, padr
from ._rmsGain import rmsGain
from ..relab import db

__all__ = ["ds_synNTFobj1", "ds_f1f2", "ds_optzeros"]


def ds_synNTFobj1(x, p, osr, f0):
    """
    Objective function for synthesizeNTF.
    """
    z = np.exp(2j*np.pi*(f0+0.5/osr*x))
    if f0 > 0:
        z = padl(z, len(p)/2, np.exp(2j*np.pi*f0))
    z = np.concatenate((z, z.conj()))
    if f0 == 0:
        z = padr(z, len(p), 1)
    f1, f2 = ds_f1f2(osr, f0)

    return db(rmsGain((z, p, 1), f1, f2))


def ds_f1f2(osr=64, f0=0, complex_flag=False):
    """
    Lower and higher extremes of the signal band as normalized frequencies

    Helper function.


    Parameters
    ----------
    osr : float, optional
        the oversamping ratio
    f0 : float, optional
        normalized center frequency for BP modulators, or 0 for LP modulators.
        Defaults to 0.
    complex_flag : bool, optional
        flag indicating if the modulator is quadrature type.

    Returns
    -------
    f1f2 : tuple with two entries corresponding to the lower and
        higher extremes of the signal band.
    """
    if complex_flag:
        f1 = f0-0.5/osr
        f2 = f0+0.5/osr
    else:
        if f0 > 0.25/osr:
            f1 = f0-0.25/osr
            f2 = f0+0.25/osr
        else:
            f1 = 0
            f2 = 0.5/osr
    return f1, f2


def ds_optzeros(n, opt=1):
    """
    Zeros which minimize the in-band noise of a delta-sigma modulator

    Helper function for synthesizeNTF, that returns the zeros as normalized
    angular frequencies


    Parameters
    ----------
    n : int
        the number of optimized zeros to return
    opt : int, optional
        flag for optimized zeros, defaults to 1
        0 -> not optimized
        1 -> optimized
        2 -> optimized with at least one zero at band-center

    Returns
    -------
    zeros : ndarray of reals
        the zeros for the modulator as normalized angular frequencies.

    Notes
    -----
    The zeros are always located on the complex unit circle. As such,
    they are returned as frequencies, not as complex values.

    The zero's frequencies are normalized with respect to the signal
    bandwidth. See also Sec. 4.3.1 in [1]_

    .. [1] Richard Schreier, Gabor C. Temes, "Understanding Delta-Sigma Data
       Converters," IEEE Press and Wiley Interscience, 2005.
    """
    if opt == 0:
        optZeros = np.zeros(ceil(n/2.))
    else:
        if n == 1:
            optZeros = np.asarray([0.])
        elif n == 2:
            if opt == 1:
                optZeros = np.asarray([sqrt(1./3)])
            else:
                optZeros = np.asarray([0.])
        elif n == 3:
            optZeros = np.asarray([sqrt(3./5), 0.])
        elif n == 4:
            if opt == 1:
                discr = sqrt(9./49-3./35)
                tmp = 3./7
                optZeros = np.sqrt([tmp+discr, tmp-discr])
            else:
                optZeros = np.asarray([0., sqrt(5./7)])
        elif n == 5:
            discr = sqrt(25./81-5./21)
            tmp = 5./9
            optZeros = np.sqrt([tmp+discr, tmp-discr, 0.])
        elif n == 6:
            if opt == 1:
                optZeros = np.asarray([0.23862059, 0.66120988, 0.9324696])
            else:
                discr = sqrt(56.)/33
                tmp = 7./11
                optZeros = np.sqrt([0, tmp+discr, tmp-discr])
        elif n == 7:
            optZeros = np.asarray([0, 0.40584371, 0.74153078, 0.94910785])
        elif n == 8:
            if opt == 1:
                optZeros = np.asarray([0.18343709, 0.52553345, 0.79666684,
                                       0.96028993])
            else:
                optZeros = np.asarray([0, 0.50563161, 0.79017286, 0.95914731])
        elif n == 9:
            optZeros = np.asarray([0, 0.32425101, 0.61337056, 0.83603082,
                                   0.9681602])
        elif n == 10:
            if opt == 1:
                optZeros = np.asarray([0.1834370913, 0.5255334458,
                                       0.7966668433, 0.9602899327])
            else:
                optZeros = np.asarray([0, 0.41572267, 0.67208682, 0.86238894,
                                       0.97342121])
        elif n == 11:
            optZeros = np.asarray([0, 0.26953955, 0.51909468, 0.73015137,
                                   0.88706238, 0.97822864])
        elif n == 12:
            if opt == 1:
                optZeros = np.asarray([0.12523875, 0.36783403, 0.58731921,
                                       0.7699033, 0.90411753, 0.9815607])
            else:
                optZeros = np.asarray([0, 0.35222363, 0.58006251, 0.76647993,
                                       0.90281326, 0.98132047])
        elif n == 13:
            optZeros = np.asarray([0, 0.23045331, 0.44849063, 0.64234828,
                                   0.8015776, 0.91759824, 0.98418306])
        elif n == 14:
            if opt == 1:
                optZeros = np.asarray([0.10806212, 0.31911586, 0.51525046,
                                       0.68729392, 0.82720185, 0.92843513,
                                       0.98628389])
            else:
                optZeros = np.asarray([0, 0.30524384, 0.50836649, 0.6836066,
                                       0.82537239, 0.92772336, 0.98615167])
        else:
            raise ValueError('Optimized zeros for n>14 are not available.')

    # Sort the zeros and replicate them.
    z = np.sort(optZeros)
    optZeros = np.zeros(n)
    m = 0
    if (n % 2) == 1:
        optZeros[0] = z[0]
        z = z[1:]
        m = m+1
    for i in xrange(len(z)):
        optZeros[m] = z[i]
        optZeros[m+1] = -z[i]
        m = m+2

    return optZeros
