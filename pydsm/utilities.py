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
Utility functions for PyDSM
===========================
"""

from __future__ import division, print_function

import numpy as np
from warnings import warn
from .exceptions import PyDsmDeprecationWarning

__all__ = ["is_negligible", "chop", "db", "cplxpair", "mdot", "EPS"]


EPS = np.finfo(float).eps


def is_negligible(x, tol=100*EPS):
    """
    Checks if a number is close to zero.

    Parameters
    ----------
    x : float, complex or array_like
        number to be checked
    tol : float, optional
        absolute tolerance. Defaults to 100 times the system epsilon.

    Returns
    -------
    y : bool or array_like of bools
        whether the input number is really close to zero or not.
    """
    return abs(np.asarray(x)) < tol


def chop(x, tol=100*EPS):
    """
    Chop to zero input numbers that are close to zero.

    Parameters
    ----------
    x : float, complex or array_like
        number to process
    tol : float, optional
        absolute tolerance. Defaults to 100 times the system epsilon.

    Returns
    -------
    y : float
        y is zero if x is close to zero according to the tolerance.
        Alternatively, y is x.

    Notes
    -----
    If the input is an array, it is always copied.
    See also `is_negligible`.
    """
    x = np.asarray(x)
    if np.iscomplexobj(x):
        return chop(x.real, tol) + 1j*chop(x.imag, tol)
    y = np.copy(x)
    y[is_negligible(y, tol)] = 0.
    return y


def mdot(*args):
    """
    Dot product taking multiple arguments

    Parameters
    ----------
    x1, ..., xn : array_like

    Returns
    -------
    y : ndarray
        y=np.dot(np.dot(np.dot(...np.dot(x1,x2),xn_2)xn_1)xn)
    """
    return reduce(np.dot, args)


def db(x, signal_type='voltage', R=1):
    """
    Converts a value to dB a la Matlab

    This function tries to replicate the ``dB`` interface of Matlab.

    .. deprecated:: 0.11.0
        Function ``db`` moved to ``relab`` module

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
    warn('Function db moved to relab module', PyDsmDeprecationWarning)
    if signal_type == 'power':
        return 10*np.log10(x)
    else:
        return 10*np.log10(np.abs(x)**2./R)


def cplxpair(x, tol=100*EPS):
    """
    Sorts values in input list by complex pairs.

    This function tries to replicate the ``cplxpair`` of Matlab, but
    currently does a very poor job.

    .. deprecated:: 0.11.0
        Function ``cplxpair`` moved to ``relab`` module

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
    warn('Function cplxpair moved to relab module',
         PyDsmDeprecationWarning)
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
