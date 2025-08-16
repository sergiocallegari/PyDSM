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
Correlation utilities (:mod:`pydsm.correlations`)
=================================================

Functions to compute the auto- and cross- correlation between two vectors.

So far only unnormalized (raw) correlations are supported.

.. currentmodule:: pydsm.correlations

Functions
---------

.. autosummary::
   :toctree: generated/

    raw_acorr  -- raw autocorrelation of a vector
    raw_xcorr  -- raw crosscorrelation of a vector
"""

from __future__ import division, print_function

import numpy as np

import sys
if sys.version_info < (3,):
    range = xrange

__all__ = ["raw_acorr", "raw_xcorr"]


def raw_acorr(x, N):
    """
    Computes the raw autocorrelation of a vector up to lag N.

    Parameters
    ----------
    x : array_like
        1-D sequence to compute the auto-correlation upon
    N : int
        the maximum (positive) lag of the raw auto-correlation to return.

    Returns
    -------
    q : ndarray
        the raw (unnormalized) autocorrelation vector.
        Assuming that :math:`m` is the length of x,
        :math:`q(k) = \\sum_{n=k}^{m-1} x(n) x(n-k)
        \\text{ for } k = 0 \\dots N`

    Notes
    -----
    The routine does not make any check on the length of x and N. It
    is responsibility of the user to assure that len(x)>=N. In some cases
    (but only in some cases), zero padding is practiced.
    """
    m = len(x)
    q = np.asarray([np.dot(x[k:m], x[0:m-k]) for k in range(N+1)])
    return q


def raw_xcorr(x, y, N):
    """
    Computes the raw crosscorrelation between two vectors up to lag N.

    Parameters
    ----------
    x : array_like
        first 1-D vector
    y : array_like
        second 1-D vector
    N : int
        the maximum (positive) lag of the raw cross-correlation to return.

    Returns
    -------
    q : ndarray
        the raw (unnormalized) crosscorrelation vector.
        Assuming that :math:`m_x` and :math:`m_y` are the lengths of x and y,
        :math:`q(k) = \\sum_{n=k}^{\\min(m_x-1,m_y+k-1)}
        x(n) y(n-k) \\text{ for } k = 0 \\dots N`

    Notes
    -----
    the routine does not make any check on the lengths of x, y and N. It
    is responsibility of the user to assure that N<=len(y). In some cases
    (but only in some cases), zero padding is assumed.
    """
    mx = len(x)
    my = len(y)
    q = np.asarray([np.dot(y[k:min(my, mx+k)],
                           x[0:min(my-k, mx)]) for k in range(N+1)])
    return q
