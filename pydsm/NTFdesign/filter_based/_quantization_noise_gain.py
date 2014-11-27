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
Computation of the noise power gain through the NTF and the output filter
=========================================================================
"""

from __future__ import division, print_function

import numpy as np
import scipy as sp
__import__("scipy.signal")
from ...ir import impulse_response
from ...delsig import evalTF
from ..weighting import quantization_weighted_noise_gain

__all__ = ["quantization_noise_gain", "quantization_noise_gain_by_conv"]


def quantization_noise_gain(NTF, H, H_type='zpk', **options):
    r"""Compute the quantization noise power gain after a filter

    Parameters
    ----------
    NTF : tuple
        NTF definition in zpk or nd form
    H : tuple or callable or array_like
        output filter definition in zpk or ba form if H_type='zpk' or 'ba'
        (in this case, H is a tuple with 3 or 2 entries);
        output filter magnitude response if H_type='mag' (in this case, H is
        a callable with argument f in [0,1/2]);
        output filter impulse response if H_type='imp' (in this case, H is an
        array)
    H_type : str
        type of specification for parameter H. One of: 'zpk', 'ba', 'mag' or
        'imp'

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

    Raises
    ------
    ValueError
        'Incorrect filter type specification' if the filter is not specified
        correctly.

    Notes
    -----
    In the default case the computation is practiced as

    .. math::
        2\int_{f=0}^{\frac{1}{2}}
        \left|H\left(\mathrm{e}^{\mathrm{i} 2\pi f}\right)\right|^2
        \left|\mathit{NTF}
        \left(\mathrm{e}^{\mathrm{i} 2\pi f}\right)\right|^2 df

     Since this function internally uses ``quantization_weighted_noise_gain``
     the latter default parameters may also affect its behavior.

    See Also
    --------
    scipy.integrate.quad : integrator used internally.
        For the meaning of the integrator parameters.
    """
    # Manage optional parameters
    opts = quantization_noise_gain.default_options.copy()
    opts.update(options)
    if H_type == 'zpk' or H_type == 'ba':
        w = lambda f: np.abs(evalTF(H, np.exp(2j*np.pi*f)))**2
    elif H_type == 'imp':
        w = lambda f: np.abs(evalTF((H, [1]), np.exp(2j*np.pi*f)))**2
    elif H_type == 'mag':
        w = lambda f: H(f)**2
    else:
        raise ValueError("Incorrect filter type specification")
    return quantization_weighted_noise_gain(NTF, w, **opts)

quantization_noise_gain.default_options = {}


def quantization_noise_gain_by_conv(NTF, H, H_type='zpk', db=80):
    """
    Computes the quantization noise power gain, based on a convolution

    Parameters
    ----------
    NTF : tuple
        NTF definition in zpk or nd form
    H : tuple or array_like
        output filter definition in zpk or nd form if H_type='zpk' or 'ba'
        (in this case, H is a tuple with 3 or 2 entries);
        output filter impulse response if H_type='imp' (in this case, H is an
        array)
    H_type : str
        type of specification for parameter H. One of: 'zpk', 'ba', or
        'imp'
    db : real
        a precision hint for the computation of impulse responses

    Returns
    -------
    a : real
        noise power gain

    Notes
    -----
    The computation is practiced as the sum of the squared entries
    in the impulse response of the cascaded filter NTF*H
    """
    h1_ir=impulse_response(NTF, db=db)
    if H_type=='zpk' or H_type=='ba':
        h2_ir=impulse_response(H, db=db)
    elif H_type=='imp':
        h2_ir = H
    else:
        raise ValueError("Incorrect filter type specification")
    conv = sp.signal.convolve(h1_ir, h2_ir)
    return np.sum(conv**2)
