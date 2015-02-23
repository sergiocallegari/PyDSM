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
Compatibility functions for PyDSM (:mod:`pydsm.relab`)
======================================================

This module re-implements some Matlab interfaces that are useful
for PyDSM.

.. currentmodule:: pydsm.relab


Functions
---------

.. autosummary::
   :toctree: generated/

   eps       -- Floating point relative accuracy
   db        -- Converts a value to dB a la Matlab
   cplxpair  -- Sorts values in input list by complex pairs
   shiftdim  -- Shift dimensions a la Matlab
"""

from __future__ import division, print_function

import numpy as np

__all__ = ["eps", "db", "cplxpair", "shiftdim"]


def eps(x):
    """Provide floating point relative accuracy

    This function tries to replicate the ``eps`` interface of Matlab.

    Parameters
    ----------
    x : float or array_like or numeric type

    Returns
    -------
    d : float or array_like
        If x is a float, the positive distance from abs(x) to the next larger
        in magnitude floating point number of the same precision as x.

        If x is an array, operation is elementwise.

        If x is a numeric type, operation is for 1.0 in that numeric type.
    """

    def _eps(xi):
        return np.finfo(xi).eps * np.abs(xi)
    if isinstance(x, (type, np.dtype)):
        return np.finfo(x).eps
    elif np.isscalar(x):
        return _eps(x)
    x = np.asarray(x)
    d = np.empty_like(x)
    d.flat = [_eps(xi) for xi in x.flat]
    return d


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


def shiftdim(x, n=None, nargout=2):
    """
    Shift dimensions a la Matlab

    When n is provided, shiftdim shifts the axes of x by n.
    If n is positive, it shifts the axes to the left, wrapping the
    leading axes with non unitary length to the end.
    When n is negative, it shifts the axes to the right, inserting n leading
    axes with unitary length.
    When n is not provided or None, it shifts the axes to the left, reducing
    the number of dimensions and removing all the leading axes with unitary
    length.

    Parameters
    ----------
    x : array like
        multi-dimensional array to operate upon
    n : int or None, optional
        amount to shift. Defaults to None, which means automatic computation
    nargout : int
        number of output values

    Returns
    -------
    y : ndarray
        the result of the axes shift operation
    n : int
        the actual shift

    Examples
    --------

    >>> from numpy.random import rand
    >>> a = rand(1, 1, 3, 1, 2)
    >>> b, n = shiftdim(a)
    >>> b.shape
    (3, 1, 2)
    >>> n
    2
    >>> c = shiftdim(b, -n, nargout=1)
    >>> np.alltrue(c == a)
    True
    >>> d = shiftdim(a, 3, nargout=1)
    >>> d.shape
    (1, 2, 1, 1, 3)

    >>> b, n = shiftdim([[[1]]])
    >>> b, n
    (array([[[1]]]), 0)
    """
    outsel = slice(nargout) if nargout > 1 else 0
    x = np.asarray(x)
    s = x.shape
    m = next((i for i, v in enumerate(s) if v > 1), 0)
    if n is None:
        n = m
    if n > 0:
        n = n % x.ndim
    if n > 0:
        if n <= m:
            x = x.reshape(s[n:])
        else:
            x = x.transpose(np.roll(range(x.ndim), -n))
    elif n < 0:
            x = x.reshape((1,)*(-n)+x.shape)
    return (x, n)[outsel]


def cplxpair(x, tol=None):
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
    x = np.asarray(x)
    if x.size == 0:
        return x
    if tol is None:
        tol = 100*eps(x.dtype)
    x = np.sort_complex(x)
    real_mask = np.abs(x.imag) < tol*np.abs(x)
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
