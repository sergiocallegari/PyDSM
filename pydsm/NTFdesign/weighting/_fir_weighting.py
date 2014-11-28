# -*- coding: utf-8 -*-

# Copyright (c) 2013, Sergio Callegari
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

These routines return the matrix from a noise weigthing function.

The matrix is symmetric Toeplitz, thus actually described by its first
row only, which is what the routines here return.
"""

from __future__ import division, print_function

import numpy as np
from ...ft import idtft_hermitian
from ...delsig import evalTF
from warnings import warn
from ...exceptions import PyDsmDeprecationWarning

__all__ = ["q0_from_noise_weighting", "q0_weighting"]


def q0_weighting(P, w, **options):
    """Compute Q matrix from a noise weighting function

    Parameters
    ----------
    P : int
        order of the FIR to be eventually synthesized
    w : callable with argument f in [0,1/2] or None or tuple
            * if function: noise weighting function
            * if None: no weighting is applied
            * if filter definition as zpk or ba tuple: weighting is implicitly
              provided by the filter

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

    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.

    See Also
    --------
    scipy.integrate.quad : integrator used internally.
        For the meaning of the integrator parameters.
    """
    # Manage parameters
    if w is None:
        w = lambda f: 1.
    elif type(w) is tuple and 2 <= len(w) <= 3:
        w = lambda f: evalTF(w, np.exp(2j*np.pi*f))**2
    # Manage optional parameters
    opts = q0_weighting.default_options.copy()
    opts.update(options)
    # Do the computation
    ac = lambda t: idtft_hermitian(w, t, **opts)
    return np.asarray(map(ac, np.arange(P+1)))

q0_weighting.default_options = {'quad_epsabs': 1E-14,
                                'quad_epsrel': 1E-9,
                                'quad_limit': 100,
                                'quad_points': None}


# Following part is deprecated

def q0_from_noise_weighting(P, w, **options):
    warn("Function superseded by q0_weighting in "
         "NTFdesign.weighting module", PyDsmDeprecationWarning)
    return q0_weighting(P, w, **options)

q0_from_noise_weighting.__doc__ = q0_weighting.__doc__ + """
    .. deprecated:: 0.11.0
    Function has been moved to the ``NTFdesign.weighting`` module with name
    ``q0_weighting``.
    """

q0_from_noise_weighting.default_options = q0_weighting.default_options
