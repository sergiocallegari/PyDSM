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

u"""
Computation of NTF related merit factors (:mod:`pydsm.NTFdesign.merit_factors`)
===============================================================================

This module provides functions for the computation of NTF merit factors

.. currentmodule:: pydsm.NTFdesign.merit_factors


Functions
---------

.. autosummary::
   :toctree: generated/

   quantization_noise_gain  -- NTF quantization noise power gain
"""

from __future__ import division, print_function

import numpy as np
from scipy.integrate import quad
from ..delsig import evalTF
from ..utilities import digested_options

__all__ = ["quantization_noise_gain"]


def quantization_noise_gain(NTF, w=None, bounds=(0, 0.5), avg=False,
                            **options):
    r"""
    Compute the NTF quantization noise power gain.

    Parameters
    ----------
    NTF : tuple
        NTF definition in zpk or nd form
    w : callable with argument f in [0,1/2] or None or tuple
            * if function: noise weighting function
            * if None: no weighting is applied
            * if filter definition as zpk or ba tuple: weighting is implicitly
              provided by the filter
    bounds : 2 elements tuple, optional
        the frequency range where the noise gain is computed. Defaults to
        (0, 0.5)
    avg: bool, optional
        If True, rather than returning the overall noise gain, the function
        returns the average noise gain over the bandwidth.

    Returns
    -------
    a : real
        noise power gain

    Other parameters
    ----------------
    quad_opts : dictionary, optional
        Parameters to be passed to the ``quad`` function used internally as
        an integrator. Allowed options are ``epsabs``, ``epsrel``, ``limit``,
        ``points``. Do not use other options since they could break the
        integrator in unexpected ways. Defaults can be set by changing the
        function ``default_options`` attribute.

    Notes
    -----
    The computation is practiced as

    .. math::
        2\int_{f=0}^{\frac{1}{2}}
        \left|\mathit{NTF}
        \left(\mathrm{e}^{\mathrm{i} 2\pi f}\right)\right|^2 w(f) df

    Use an on-off weighting function :math:`w(f)` for multiband evaluation.

    In case the weighting function has discontinuities, report them to the
    integrator via the ``quad_points`` parameter.

    See Also
    --------
    scipy.integrate.quad : for the meaning of the integrator parameters.
    """
    # Manage parameters
    if w is None:
        w = lambda f: 1.
    elif type(w) is tuple and 2 <= len(w) <= 3:
        h = w
        w = lambda f: np.abs(evalTF(h, np.exp(2j*np.pi*f)))**2
    # Manage optional parameters
    opts = digested_options(options, quantization_noise_gain.default_options,
                            [], ['quad_opts'])
    # Compute
    c = 1/(bounds[1]-bounds[0]) if avg else 2.
    return c*quad(lambda f: np.abs(evalTF(NTF, np.exp(2j*np.pi*f)))**2*w(f),
                  bounds[0], bounds[1], **opts["quad_opts"])[0]

quantization_noise_gain.default_options = {"quad_opts": {"epsabs": 1E-14,
                                                         "epsrel": 1E-9,
                                                         "limit": 100,
                                                         "points": None}}
