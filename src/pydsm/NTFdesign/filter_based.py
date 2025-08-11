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
Output filter based NTF synthesis (:mod:`pydsm.NTFDesign.filter_based`)
=======================================================================

This modules provides code for the synthesis of the modulator NTF, based
on the filter placed after the modulator for the removal of the quantization
noise.

.. deprecated:: 0.11.0
    Use `ntf_fir_weighting` or functions from `NTFdesign.weighting`
    module.

.. currentmodule:: pydsm.NTFdesign.filter_based


Deprecated functions
--------------------

.. autosummary::
   :toctree: generated/

   quantization_noise_gain -- Alternate quantization noise power gain func
   quantization_noise_gain_by_conv -- Alternate function based on convolution
   q0_from_filter -- Alternate compute Q matrix from the output filter
   synthesize_ntf_from_filter -- Alternate compute FIR NTF from output filter
"""

from __future__ import division, print_function

from .merit_factors import quantization_noise_gain as _quantization_noise_gain
from .legacy import (quantization_noise_gain_by_conv as
                     _quantization_noise_gain_by_conv,
                     q0_from_filter_ir)
from .weighting import q0_weighting, ntf_fir_from_q0
from warnings import warn
from ..exceptions import PyDsmDeprecationWarning
from ..utilities import digested_options

__all__ = ["quantization_noise_gain", "quantization_noise_gain_by_conv",
           "q0_from_filter", "synthesize_ntf_from_filter"]


def quantization_noise_gain(NTF, H, H_type='zpk', **options):
    r"""Compute the quantization noise power gain after a filter

    .. deprecated:: 0.11.0
       Use :func:`pydsm.NTFdesign.quantization_noise_gain`. Note that the
       interface is slightly different.

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
    quad_opts : dictionary, optional
        Parameters to be passed to the ``quad`` function used internally as
        an integrator. Allowed options are ``epsabs``, ``epsrel``, ``limit``,
        ``points``. Do not use other options since they could break the
        integrator in unexpected ways. Defaults can be set by changing the
        function ``default_options`` attribute.

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

    See Also
    --------
    scipy.integrate.quad : for the meaning of the integrator parameters.
    """
    warn("Function meant for removal", PyDsmDeprecationWarning)
    # Manage optional parameters
    opts = digested_options(options, quantization_noise_gain.default_options,
                            [], ['quad_opts'])
    if H_type == 'zpk' or H_type == 'ba':
        w = H
    elif H_type == 'imp':
        w = (H, [1])
    elif H_type == 'mag':
        w = lambda f: H(f)**2
    else:
        raise ValueError("Incorrect filter type specification")
    return _quantization_noise_gain(NTF, w, **opts)

quantization_noise_gain.default_options = \
    _quantization_noise_gain.default_options.copy()


def quantization_noise_gain_by_conv(NTF, H, H_type='zpk', db=80):
    """
    Alias of :func:`pydsm.NTFdesign.legacy.quantization_noise_gain_by_conv`

    .. deprecated:: 0.11.0
        Function has been moved to the :mod"`NTFdesign.legacy` module.
    """
    warn("Function moved in ``NTFdesign.legacy`` module",
         PyDsmDeprecationWarning)
    return _quantization_noise_gain_by_conv(NTF, H, H_type, db)


def q0_from_filter(P, H, H_type='zpk', **options):
    """
    Compute Q matrix from the modulator output filter

    .. deprecated:: 0.11.0
       Use :func:`NTFdesign.weighting.q0_weighting`. Note that the
       interface is slightly different.

    Parameters
    ----------
    P : int
        order of the FIR to be eventually synthesized
    H : tuple or array_like or callable
        output filter description.
        This is given by a zpk or ba form if F_type is 'zpk' or 'ba'.
        It is a magnitude response (function of f, with f in [0,1/2]) if
        F_type is 'mag'.
        It is an impulse response if F_type is 'imp'.
    H_type : str
        string indicating how F is expressed. Can be 'zpk', 'ba', 'mag' or
        'imp'

    Returns
    -------
    q0 : ndarray
        the first row of the matrix Q used in the NTF optimization

    Other parameters
    ----------------
    quad_opts : dictionary, optional
        Parameters to be passed to the ``quad`` function used internally as
        an integrator. Allowed options are ``epsabs``, ``epsrel``, ``limit``,
        ``points``. Do not use other options since they could break the
        integrator in unexpected ways. Defaults can be set by changing the
        function ``default_options`` attribute.

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
    scipy.integrate.quad : for the meaning of the integrator parameters.
    """
    warn("Function meant for removal", PyDsmDeprecationWarning)
    # Manage optional parameters
    opts = digested_options(options, q0_from_filter.default_options,
                            [], ['quad_opts'])
    # Do the computation
    if H_type == 'zpk' or H_type == 'ba':
        w = H
        q0 = q0_weighting(P, w, **opts)
    elif H_type == 'imp':
        q0 = q0_from_filter_ir(P, H)
    elif H_type == 'mag':
        w = lambda f: H(f)**2
        q0 = q0_weighting(P, w, **opts)
    else:
        raise ValueError("Incorrect filter type specification")
    return q0

q0_from_filter.default_options = q0_weighting.default_options.copy()


def synthesize_ntf_from_filter(order, H, H_type='zpk', H_inf=1.5,
                               normalize="auto", **options):
    u"""Synthesize a FIR NTF based on the ΔΣM output filter.

    The ΔΣ modulator NTF is designed after a specification of the
    filter in charge of removing the quantization noise

    .. deprecated:: 0.11.0
       Use :func:`NTFdesign.ntf_fir_weighting`. Note that the interface
       is slightly different.

    Parameters
    ----------
    order : int
        Delta sigma modulator order
    F : tuple or callable or array_like
        filter specification, the format depends on parameter F_type.
        a zpk or ba tuple if F_type is 'zpk' or 'ba', respectively.
        a function of f, for f in [0,1/2] if F_type is 'mag'
        an array containing an impulse response if F_type is 'imp'
    F_type : str
        string indicating the type of filter specification. Can be 'zpk',
        'ba', 'mag' or 'imp'.
    H_inf : real, optional
        Max peak NTF gain, defaults to 1.5, used to enforce the Lee criterion
    normalize : string or real, optional
        Normalization to apply to the quadratic form used in the NTF
        selection. Defaults to 'auto' which means setting the top left entry
        in the matrix Q defining the quadratic form to 1.

    Returns
    -------
    ntf : ndarray
        FIR NTF in zpk form

    Other parameters
    ----------------
    show_progress : bool, optional
        provide extended output, default is True and can be updated by
        changing the function ``default_options`` attribute.
    cvxopt_opts : dictionary, optional
        A dictionary of options for the ``cvxopt`` optimizer
        Allowed options include:

        ``maxiters``
            Maximum number of iterations (defaults to 100)
        ``abstol``
            Absolute accuracy (defaults to 1e-7)
        ``reltol``
            Relative accuracy (defaults to 1e-6)
        ``feastol``
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break `cvxopt` in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.
    quad_opts : dictionary, optional
        Parameters to be passed to the ``quad`` function used internally as
        an integrator. Allowed options are ``epsabs``, ``epsrel``, ``limit``,
        ``points``. Do not use other options since they could break the
        integrator in unexpected ways. Defaults can be set by changing the
        function ``default_options`` attribute.

    See Also
    --------
    scipy.integrate.quad : for the meaning of the integrator parameters.
    cvxopt : for the optimizer parameters
    """
    warn("Function meant for removal", PyDsmDeprecationWarning)
    # Manage optional parameters
    opts1 = digested_options(options,
                             synthesize_ntf_from_filter.default_options,
                             ['quad_opts'], [], False)
    opts2 = digested_options(options,
                             synthesize_ntf_from_filter.default_options,
                             ['show_progress'], ['cvxopt_opts'])
    # Do the computation
    q0 = q0_from_filter(order, H, H_type, **opts1)
    return ntf_fir_from_q0(q0, H_inf, normalize, **opts2)

synthesize_ntf_from_filter.default_options = \
    q0_from_filter.default_options.copy()
synthesize_ntf_from_filter.default_options.update(
    ntf_fir_from_q0.default_options)
