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
Computation of the impulse response of a DT filter (:mod:`pydsm.ir`)
====================================================================

Compute (approximating by truncation) the impulse response
of a discrete time filter, trying to (roughly) guess the appropriate
truncation.

.. currentmodule:: pydsm.ir

Functions
---------

.. autosummary::
   :toctree: generated/

    guess_ir_length  -- Guess appropriate truncation length
    impulse_response  -- Compute impulse response of DT filter
"""

from __future__ import division, print_function

import numpy as np
import scipy as sp
__import__("scipy.signal")

__all__ = ["impulse_response", "guess_ir_length"]


def impulse_response(h, m=None, db=80):
    """
    Computes the filter impulse response

    Parameters
    ----------
    h : tuple_like
        the filter definition either in zpk or in nd form.

    Returns
    -------
    ir : ndarray
        the truncated impulse response

    Other Parameters
    ----------------
    m : int, optional
        the number of samples after which the impulse response should be
        truncated. Defaults to None, which means *try to guess*
    db : real, optional
        a hint about how to guess the length where the impuls response
        should be truncated (defaults to 80)

    Notes
    -----
    The guess about the lenght where the impulse response can be truncated
    is extremely rough. See :func:`guess_ir_length` in this module for
    further info.
    """
    if len(h) == 3:
        (b, a) = sp.signal.zpk2tf(*h)
        b = b.real
        a = a.real
    else:
        (b, a) = h
    if m is None:
        m = guess_ir_length(h, db)
    ins = np.zeros(m)
    ins[0] = 1
    return sp.signal.lfilter(b, a, ins)


def guess_ir_length(h, db=80):
    """
    Tries to estimate an appropriate length for the filter response

    Parameters
    ----------
    h : tuple_like
        the filter definition either in zpk or in nd form.
    db : real, optional
        a hint about how to guess the length where the impulse response
        should be truncated. This is defined on a log scale. The larger
        the longer the resulting length. Defaults to 80.

    Returns
    -------
    m : int
        a guess about the appropriate number of samples to represent the
        filter impulse response with the required accuracy

    Notes
    -----
    The guess is based on the slowlest pole of the filter, considering when
    its response is attenuated to -db. This can be by far too optimistic in
    some cases and particularly when there are overlapping or similar poles.

    Do not try to use this function for filters with poles in 1.
    """
    # Put h in zpk form if it is in tf form
    if len(h) == 2:
        h = sp.signal.tf2zpk(*h)
    pp = h[1]
    t_z = len(h[0])+1
    if len(pp) == 0:
        t_p = 0
    else:
        # Try to estimate length of decay of the filter h.
        # The estimation is extremely rough, based on the decay
        # rate of the pole with maximum magnitude.
        # Thus, it breaks easily when there are poles very close
        # one to the other or overlapping.
        # Furthermore, this code should not be called if there is
        # a pole in z=1.
        os = np.seterr(divide='ignore')
        sr = np.log(np.abs(pp))
        np.seterr(**os)
        # Take slowlest pole
        wmin = np.min(np.abs(sr))
        # 1/omega min is time constant in sample periods.
        # Let's multiply the time constant in order to have
        # the transient attenuated by db decibels
        t_p = int(np.ceil(db/20*np.log(10)/wmin))
    return t_p+t_z
