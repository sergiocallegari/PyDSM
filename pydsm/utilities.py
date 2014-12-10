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
from .import relab

__all__ = ["is_negligible", "chop", "db", "cplxpair", "mdot", "EPS",
           "split_options", "strip_options"]


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


def split_options(opts, prefixes=[], keys=[]):
    """
    Splits a set of keyword arguments by name and prefix.
    """
    o = {}
    for prefix in prefixes:
        o[prefix] = {}
    for key in opts:
        found = False
        if key in keys:
            o[key] = opts[key]
        else:
            for prefix in prefixes:
                if key.startswith(prefix):
                    o[prefix][key] = opts[key]
                    found = True
                    break
            if not found:
                raise TypeError("Unexpected keyword argument '%s'" % key)
    return o


def strip_options(opts, prefix):
    """
    Strips a prefix from a set of keyword arguments.
    """
    o = {}
    opts = opts[prefix]
    for key in opts:
        o[key[len(prefix):]] = opts[key]
    return o


# Following two functions are deprecated

def db(x, signal_type='voltage', R=1):
    warn('Function db moved to relab module', PyDsmDeprecationWarning)
    return relab.db(x, signal_type, R)

db.__doc__ = relab.db.__doc__ + """
     .. deprecated:: 0.11.0
        Function ``db`` moved to ``relab`` module
"""


def cplxpair(x, tol=100*EPS):
    warn('Function cplxpair moved to relab module',
         PyDsmDeprecationWarning)
    return relab.cplxpair(x, tol)

cplxpair.__doc__ = relab.cplxpair.__doc__ + """
    .. deprecated:: 0.11.0
        Function ``cplxpair`` moved to ``relab`` module
"""
