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
Collection of DELSIG style decibel routines
===========================================
"""

import numpy as np

__all__=["dbv", "dbp", "dbm", "undbv", "undbp", "undbm"]

def dbv(x):
    """
    Converts a voltage (or current) ratio to dB.

    Parameters
    ----------
    x : real or array_like of reals
        argument should be strictly positive

    Returns
    -------
    y : real or array_like of reals
        20*log10(x)
    """
    return 20*np.log10(x)

def dbp(x):
    """
    Converts a power ratio to dB.

    Parameters
    ----------
    x : real or array_like of reals
        argument should be strictly positive

    Returns
    -------
    y : real or array_like of reals
        10*log10(x)
    """
    return 10*np.log10(x)

def dbm(v, R=50.):
    """
    Converts argument from rms voltage to dBm.

    Parameters
    ----------
    v : real or array_like of reals
        rms value to be converted in dBm
    R : real, optional
        value of test resistor (defaults to 50)

    Returns
    -------
    y : real or array_like of reals
        dBm value corresponding to the power provided by the rms voltag v
        on the test resistor R, normalized over the reference value of 1mW
    """
    return 10*np.log10(np.abs(v**2)/R)+30

def undbv(x):
    """
    Converts argument from db to a voltage ratio.

    Inverse of `dbv`.

    Parameters
    ----------
    x : real or array_like of reals
        value to convert

    Returns
    -------
    y : real or array_like of reals
        10^(x/20)
    """
    return 10**(np.asarray(x)/20.)

def undbp(x):
    """
    Converts argument from db to a power ratio.

    Inverse of `dbp`.

    Parameters
    ----------
    x : real or array_like of reals
        value to convert

    Returns
    -------
    y : real or array_like of reals
        10^(x/10)
    """
    return 10**(np.asarray(x)/10.)

def undbm(x, R=50.):
    """
    Converts argument from dBm power to rms voltage.

    Inverse of `dbm`.

    Parameters
    ----------
    x : real or array_like of reals
        value to convert
    R : real, optional
        test resistor

    Returns
    -------
    y : real or array_like of reals
        rms voltage that applied over a resistor R gives a power equal to x
        in dBm (i.e. referred to 1mW)
    """
    return np.sqrt(R*10**((x-30)/10.))
