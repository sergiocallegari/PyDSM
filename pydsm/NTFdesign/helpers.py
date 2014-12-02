# -*- coding: utf-8 -*-

# Copyright (c) 2014, Sergio Callegari
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

u"""
Helper functions for NTF design (:mod:`pydsm.NTFdesign.helpers`)
================================================================

This modules provides helper functions for the design of the NTF
of ΔΣ modulators.

.. currentmodule:: pydsm.NTFdesign.helpers

.. autosummary::
   :toctree: generated/

   maxflat_fir_zeros            -- Zeros of a maxflat FIR transfer function
"""


from __future__ import division

import numpy as np

__all__ = ["maxflat_fir_zeros"]


def maxflat_fir_zeros(order, alpha):
    """
    Compute the zeros of a maxflat FIR transfer function.

    The computed FIR transfer function is used in the DELSIG
    :func:`pydsm.delsig.synthesizeDSM()` design method
    (also known as :func:`pydsm.NTFdesign.ntf_schreier()`)
    for the denominator of the noise transfer function.

    Parameters
    ----------
    order : int
        the transfer function order
    alpha : float
        parameter controlling the steepness of the magnitude
        response

    Returns
    -------
    zeros : ndarray
        an array of complex roots of the FIR transfer function
    """
    x = 1./np.sqrt(alpha)
    me2 = -0.5*(x**(2./order))
    w = (2*np.arange(1, order+1)+1)*np.pi/order
    mb2 = 1+me2*np.exp(1j*w)
    p = mb2 - np.sqrt(mb2**2-1)
    # Reflect poles that fall out of the unit circle
    out = (np.abs(p) > 1)
    p[out] = 1/p[out]
    return p
