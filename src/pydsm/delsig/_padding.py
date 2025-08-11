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
DELSIG style matrix and vector padding routines
===============================================
"""

import numpy as np

__all__ = ["padl", "padr", "padt", "padb"]


def padl(x, n, val=0):
    """
    Pad a matrix or a vector on the left.

    Parameters
    ----------
    x : ndarray
        matrix or vector to be padded (must be 1D or 2D)
    n : int
        (column) length to pad x to
    val : scalar, optional
        fill value for padding. Defaults to 0.

    Returns
    -------
    y : ndarray
        matrix or vector x, padded to the left to have n entries (if 1D)
        or n columns (if 2D). Padding value is val.

    Notes
    -----
    If x is an array rather than a matrix, then an array is returned,
    padded at its beginning. Namely, the array is interpreted as a row
    vector.
    """
    if x.ndim == 1:
        return np.hstack((val*np.ones(n-x.shape[0]), x))
    return np.hstack((val*np.ones((x.shape[0], n-x.shape[1])), x))


def padr(x, n, val=0):
    """Pad a matrix or a vector x on the right.

    Parameters
    ----------
    x : ndarray
        matrix or vector to be padded (must be 1D or 2D)
    n : int
        (column) length to pad x to
    val : scalar, optional
        fill value for padding. Defaults to 0.

    Returns
    -------
    y : ndarray
        matrix or vector x, padded to the right to have n entries (if 1D)
        or n columns (if 2D). Padding value is val.

    Notes
    -----
    If x is an array rather than a matrix, then an array is returned,
    padded at its end. Namely, the array is interpreted as a row
    vector.
    """
    if x.ndim == 1:
        return np.hstack((x, (val*np.ones(n-x.shape[0]))))
    return np.hstack((x, (val*np.ones((x.shape[0], n-x.shape[1])))))


def padt(x, n, val=0):
    """Pad a matrix or a vector on the top.

    Parameters
    ----------
    x : ndarray
        matrix or vector to be padded (must be 1D or 2D)
    n : int
        row height to pad x to
    val : number, optional
        fill value for padding. Defaults to 0.

    Returns
    -------
    y : ndarray
        matrix or vector x, padded to the top to have n rows.
        Padding value is val.

    Notes
    -----
    If x is an array rather than a matrix, then an array is returned,
    padded at its beginning. Namely, the array is interpreted as a column
    vector.
    """
    if x.ndim == 1:
        return np.hstack((val*np.ones(n-x.shape[0]), x))
    return np.vstack((val*np.ones((n-x.shape[0], x.shape[1])), x))


def padb(x, n, val=0):
    """Pad a matrix or a vector on the bottom.

    Parameters
    ----------
    x : ndarray
        matrix or vector to be padded (must be 1D or 2D)
    n : int
        row height to pad x to
    val : scalar, optional
        fill value for padding. Defaults to 0.

    Returns
    -------
    y : ndarray
        matrix or vector x, padded to the bottom to have n rows.
        Padding value is val.

    Notes
    -----
    If x is an array rather than a matrix, then an array is returned,
    padded at its end. Namely, the array is interpreted as a column
    vector.
    """
    if x.ndim == 1:
        return np.hstack((x, (val*np.ones(n-x.shape[0]))))
    return np.vstack((x, val*np.ones((n-x.shape[0], x.shape[1]))))
