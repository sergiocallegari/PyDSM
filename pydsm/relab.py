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


def eps(x=1.):
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

    Examples
    --------
    >>> eps()
    2.2204460492503131e-16
    >>> eps(1/2)
    1.1102230246251565e-16
    >>> eps(np.float32)
    1.1920929e-07
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


def cplxpair(x, tol=None, dim=None):
    """
    Sorts values into complex pairs a la Matlab.

    The function takes a vector or multidimensional array of of complex
    conjugate pairs or real numbers and rearranges it so that the complex
    numbers are collected into matched pairs of complex conjugates. The pairs
    are ordered by increasing real part, with purely real elements placed
    after all the complex pairs.

    In the search for complex conjugate pairs a relative tolerance equal to
    ``tol`` is used for comparison purposes. The default tolerance is
    100 times the system floating point accuracy.

    If the input vector is a multidimensional array, the rearrangement is done
    working along the axis specifid by the parameter ``dim`` or along the
    first axis with non-unitary length if ``dim`` is not provided.

    Parameters
    ----------
    x : array_like of complex
        x is an array of complex values, with the assumption that it contains
        either real values or complex values in conjugate pairs.
    tol: real, optional
        relative tolerance for the recognition of pairs.
        Defaults to 100 times the system floating point accuracy for the
        specific number type.
    dim: integer, optional
        The axis to operate upon.

    Returns
    -------
    y : ndarray
        y is an array of complex values, with the same values in x, yet now
        sorted as complex pairs by increasing real part. Real elements in x
        are place after the complex pairs, sorted in increasing order.

    Raises
    ------
    ValueError
        'Complex numbers cannot be paired' if there are unpaired complex
        entries in x.

    Examples
    --------
    >>> a = np.exp(2j*np.pi*np.arange(0, 5)/5)
    >>> b1 = cplxpair(a)
    >>> b2 = np.asarray([-0.80901699-0.58778525j, -0.80901699+0.58778525j,
    ...                   0.30901699-0.95105652j,  0.30901699+0.95105652j,
    ...                   1.00000000+0.j])
    >>> np.allclose(b1, b2)
    True

    >>> cplxpair(1)
    array([1])

    >>> cplxpair([[5, 6, 4], [3, 2, 1]])
    array([[3, 2, 1],
           [5, 6, 4]])

    >>> cplxpair([[5, 6, 4], [3, 2, 1]], dim=1)
    array([[4, 5, 6],
           [1, 2, 3]])

    See also
    --------
    eps : the system floating point accuracy
    """

    def cplxpair_vec(x, tol):
        real_mask = np.abs(x.imag) <= tol*np.abs(x)
        x_real = np.sort(np.real(x[real_mask]))
        x_cplx = np.sort(x[np.logical_not(real_mask)])
        if x_cplx.size == 0:
            return x_real
        if (x_cplx.size % 2) != 0:
            raise ValueError('Complex numbers cannot be paired')
        if np.any(np.real(x_cplx[1::2])-np.real(x_cplx[0::2]) >
                  tol*np.abs(x_cplx[0::2])):
            raise ValueError('Complex numbers cannot be paired')
        start = 0
        while start < x_cplx.size:
            sim_len = next((i for i, v in enumerate(x_cplx[start+1:]) if
                           (np.abs(np.real(v)-np.real(x_cplx[start])) >
                            tol*np.abs(v))), x_cplx.size-start-1)+1
            if (sim_len % 2) != 0:
                sim_len -= 1
            # At this point, sim_len elements with identical real part
            # have been identified.
            sub_x = x_cplx[start:start+sim_len]
            srt = np.argsort(np.imag(sub_x))
            sub_x = sub_x[srt]
            if np.any(np.abs(np.imag(sub_x)+np.imag(sub_x[::-1])) >
                      tol*np.abs(sub_x)):
                raise ValueError('Complex numbers cannot be paired')
            # Output should contain "perfect" pairs. Hence, keep entries
            # with positive imaginary parts amd use conjugate for pair
            x_cplx[start:start+sim_len] = np.concatenate(
                (np.conj(sub_x[:sim_len//2-1:-1]),
                 sub_x[:sim_len//2-1:-1]))
            start += sim_len
        return np.concatenate((x_cplx, x_real))

    x = np.atleast_1d(x)
    if x.size == 0:
        return x
    if dim is None:
        dim = next((i for i, v in enumerate(x.shape) if v > 1), 0)
    if tol is None:
        try:
            tol = 100*eps(x.dtype)
        except:
            tol = 100*eps(np.float)
    return np.apply_along_axis(cplxpair_vec, dim, x, tol)
