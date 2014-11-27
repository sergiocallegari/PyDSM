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
Functions to obtain the matrix used in the NTF optimization
===========================================================

The matrix is symmetric Toeplitz, thus actually described by its first
row only, which is what the routines here return.
"""

from __future__ import division, print_function

from ...correlations import raw_acorr
from ..weighting import q0_from_noise_weighting
from ...delsig import evalTF
import numpy as np

__all__ = ["q0_from_filter"]


def q0_from_filter(P, F, F_type='zpk', **options):
    """Compute Q matrix from the modulator output filter

    Parameters
    ----------
    P : int
        order of the FIR to be eventually synthesized
    F : tuple or array_like or callable
        output filter description.
        This is given by a zpk or ba form if F_type is 'zpk' or 'ba'.
        It is a magnitude response (function of f, with f in [0,1/2]) if
        F_type is 'mag'.
        It is an impulse response if F_type is 'imp'.
    F_type : str
        string indicating how F is expressed. Can be 'zpk', 'ba', 'mag' or
        'imp'

    Returns
    -------
    q0 : ndarray
        the first row of the matrix Q used in the NTF optimization

    Other parameters
    ----------------
    quad_xxx : various type
        Parameters prefixed by ``quad_`` are passed to the ``quad``
        function that is used internally as an integrator. Allowed options
        are ``quad_epsabs``, ``quad_epsrel``, ``quad_limit``, ``quad_points``.
        Do not use other options since they could break the integrator in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    Raises
    ------
    ValueError
        'Incorrect filter type specification' if the filter is not specified
        correctly.

    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.

    See Also
    --------
    scipy.integrate.quad : integrator used internally.
        For the meaning of the integrator parameters.
    """
    # Manage optional parameters
    opts = q0_from_filter.default_options.copy()
    opts.update(options)
    # Do the computation
    if F_type == 'zpk' or F_type == 'ba':
        w = lambda f: np.abs(evalTF(F, np.exp(2j*np.pi*f)))**2
        q0 = q0_from_noise_weighting(P, w, **opts)
    elif F_type == 'imp':
        q0 = raw_acorr(F, P)
    elif F_type == 'mag':
        w = lambda f: F(f)**2
        q0 = q0_from_noise_weighting(P, w, **opts)
    else:
        raise ValueError("Incorrect filter type specification")
    return q0

q0_from_filter.default_options = \
    q0_from_noise_weighting.default_options.copy()
