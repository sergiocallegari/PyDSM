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

"""
Compatibility functions for PyDSM
=================================

This module re-implements some matlab interfaces that are useful
for PyDSM.
"""

from __future__ import division, print_function

import numpy as np
from .utilities import chop


def db(x, signal_type='voltage', R=1):
    """
    Converts a value to dB a la Matlab

    This function tries to replicate the ``dB`` interface of Matlab.

    Parameters
    ----------
    x : real
        value to be converted. Should be positive and non null.
    signal_type : string, optional
        either 'voltage' or 'power'. Defaults to 'voltage'
    R : real, optional
        load resistance value. Defaults to 1.

    Returns
    -------
    y : real
        value in dB corresponding to x
        if signal_type is 'power' the result is 10*log10(x). Otherwise (if
        signal_type is 'voltage'), then power is measured over resistor R

    Notes
    -----
    The default R value assures that when signal_type is 'voltage' dB defaults
    to the classical 20*log10(x) computation.
    """
    if signal_type == 'power':
        return 10*np.log10(x)
    else:
        return 10*np.log10(np.abs(x)**2./R)


def cplxpair(x, tol=100*np.finfo(float).eps):
    """
    Sorts values in input list by complex pairs.

    This function tries to replicate the ``cplxpair`` of Matlab, but
    currently does a very poor job.

    Parameters
    ----------
    x : array_like of complex
        x is an array of complex values, with the assumption that it contains
        either real values or complex values in conjugate couples.
    tol: real, optional
        absolute tolerance for the recognition of pairs.
        Defaults to 100 times the system epsilon.

    Returns
    -------
    y : ndarray
        y is an array of complex values, with the same values in x, yet now
        sorted as complex pairs by increasing real part. Real elements in x
        are place after the complex pairs, sorted in increasing order.

    Raises
    ------
    ValueError
        'Cannot identify complex pairs.' if there are unpaired complex entries
        in x.

    Notes
    -----
    This function is similar to Matlab cplxpair, but not quite.

    """
    x = np.sort_complex(chop(x, tol))
    real_mask = np.isreal(x)
    x_real = x[real_mask]
    x_cplx = x[np.logical_not(real_mask)]
    pos_mask = x_cplx.imag > 0
    x_cplx_pos = x_cplx[pos_mask]
    x_cplx_neg = x_cplx[np.logical_not(pos_mask)]
    x_cplx2 = 1j*np.zeros(2*x_cplx_pos.size)
    for i in range(len(x_cplx_pos)):
        x_cplx2[2*i] = x_cplx_pos[i]
        x_cplx2[2*i+1] = x_cplx_neg[i]
        if abs(x_cplx2[2*i]-np.conj(x_cplx2[2*i+1])) > tol:
            raise ValueError('Cannot identify complex pairs.')
    return np.concatenate((x_cplx2, x_real))
