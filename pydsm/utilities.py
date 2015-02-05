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
Generic utility functions (:mod:`pydsm.utilities`)
==================================================

This module includes some generic utility functions. The statbility of these
interfaces may not be guaranteed.

.. currentmodule:: pydsm.utilities


Constants
---------

.. autodata:: EPS
   :annotation:


Functions
---------

.. autosummary::
   :toctree: generated/

   is_negligible     -- Check if a number is close to zero
   chop              -- Chop to zero numbers that are close to zero
   mdot              -- Dot product taking multiple arguments
   digested_options  -- Helper function for the management of default options

Deprecated functions
--------------------

.. autosummary::
   :toctree: generated/

   db        -- Alias for :func:`relab.db`
   cplxpair  -- Alias for :func:`relab.cplxpair`
"""

from __future__ import division, print_function

import numpy as np
from warnings import warn
from .exceptions import PyDsmDeprecationWarning
from . import relab

__all__ = ["is_negligible", "chop", "db", "cplxpair", "mdot", "EPS",
           "digested_options"]


EPS = np.finfo(float).eps
"""The value returned by ``np.finfo(float).eps``, namely the smallest
   representable positive numebr such that ``1.0 + EPS != 1``."""


def is_negligible(x, tol=100*EPS):
    """
    Check if a number is close to zero.

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
    Dot product taking multiple arguments.

    Parameters
    ----------
    x1, ..., xn : array_like

    Returns
    -------
    y : ndarray
        y=np.dot(np.dot(np.dot(...np.dot(x1,x2),xn_2)xn_1)xn)

    Notes
    -----
    This function is not necessary with recent numpy version where
    one can write things like ``A.dot(B).dot(C)``
    """
    return reduce(np.dot, args)


def digested_options(opts, defaults, keys=[], multikeys=[], emptycheck=True):
    """
    Helper function for the management of default options

    This is used to validate and pre-digest the parameters of functions
    taking a variable number of keyword arguments with the ``**kwargs``
    syntax.

    Parameters
    ----------
    opts : dict
        Dictionary of optional keyword arguments
    defaults : dict
        Dictionary of defaul values
    keys : list
        List of keys to extract whose value can be directly updated from the
        defaults
    multikeys : list
        List of keys to extract whose value must be recursively updated from
        the defaults dictionary (namely, value of a multikey is a dictionary
        itself)
    emptychec : boolean
        Whether to rise an error if items remain in opts after processing.

    Returns
    -------
    eopts : dict
        Dictionary of extracted and validated options. This dictionary
        contains the values from the defaults dictionary, updated with
        the values passed in opts.

    Notes
    -----
    On exit, the input dictionary opt is changed, with all the validated
    options removed from it.
    """
    out = {}
    for key in keys:
        out[key] = opts.pop(key, defaults[key])
    for key in multikeys:
        d = defaults[key].copy()
        d.update(opts.pop(key, []))
        out[key] = d
    if emptycheck and (len(opts) != 0):
        raise TypeError("Unexpected keyword argument(s) '%s'" %
                        list(opts.viewkeys()))
    return out


# Following two functions are deprecated

def db(x, signal_type='voltage', R=1):
    """
    Alias for :func:`pydsm.relab.db`.

     .. deprecated:: 0.11.0
        Function :func:`db` moved to :mod:`pydsm.relab` module
    """
    warn('Function db moved to relab module', PyDsmDeprecationWarning)
    return relab.db(x, signal_type, R)


def cplxpair(x, tol=100*EPS):
    """
    Alias for :func:`pydsm.relab.cplxpair`.

     .. deprecated:: 0.11.0
        Function :func:`cplxpair` moved to :mod:`pydsm.relab` module
    """
    warn('Function cplxpair moved to relab module',
         PyDsmDeprecationWarning)
    return relab.cplxpair(x, tol)
