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
Computation of NTF related merit factors
========================================
"""

from __future__ import division, print_function

import numpy as np
from scipy.integrate import quad
from ..delsig import evalTF

__all__ = ["quantization_noise_gain"]


def quantization_noise_gain(NTF, w=None, bounds=(0, 0.5), **options):
    r"""Computes the NTF quantization noise power gain

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

    Returns
    -------
    a : real
        noise power gain

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
    scipy.integrate.quad : integrator used internally.
        For the meaning of the integrator parameters.
    """
    # Manage parameters
    if w is None:
        w = lambda f: 1.
    elif type(w) is tuple and 2 <= len(w) <= 3:
        w = lambda f: evalTF(w, np.exp(2j*np.pi*f))**2
    # Manage optional parameters
    opts = quantization_noise_gain.default_options.copy()
    opts.update(options)
    quad_opts = {k[5:]: v for k, v in opts.iteritems()
                 if k.startswith('quad_')}
    # Compute
    return 2*quad(lambda f: np.abs(evalTF(NTF, np.exp(2j*np.pi*f)))**2*w(f),
                  bounds[0], bounds[1], **quad_opts)[0]

quantization_noise_gain.default_options = {'quad_epsabs': 1.49e-08,
                                           'quad_epsrel': 1.49e-08,
                                           'quad_limit': 50,
                                           'quad_points': None}
