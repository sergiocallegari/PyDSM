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
Synthesize a FIR NTF from specs of filter used to remove quantization noise
===========================================================================
"""

from __future__ import division, print_function

from ._q0_from_filter import q0_from_filter
from ..weighting import synthesize_ntf_from_q0
from warnings import warn
from ...errors import PyDsmWarning

__all__ = ["synthesize_ntf_from_filter",
           "synthesize_ntf_from_filter_imp",
           "synthesize_ntf_from_filter_mag",
           "synthesize_ntf_from_filter_ir",
           "synthesize_ntf_from_q0"]


def synthesize_ntf_from_filter_ir(order, h_ir, H_inf=1.5, normalize="auto",
                                  **options):
    warn('Deprecated function synthesize_ntf_from_filter_ir.\n'
         'Will be removed shortly.\n'
         'Use synthesize_ntf_from_filter instead.', PyDsmWarning)
    # Manage optional parameters
    opts = synthesize_ntf_from_filter_ir.default_options.copy()
    opts.update(options)
    # Do the computation
    return synthesize_ntf_from_filter(order, h_ir, 'imp', H_inf, normalize,
                                      *opts)

synthesize_ntf_from_filter_ir.default_options = {'quad_epsabs': 1E-14,
                                                 'quad_epsrel': 1E-9}


def synthesize_ntf_from_filter_imp(order, h_ir, H_inf=1.5, normalize="auto",
                                   **options):
    u"""Synthesize a FIR NTF based on the ΔΣM output filter impulse response.

    The ΔΣ modulator NTF is designed after the impulse response of the filter
    in charge of removing the quantization noise

    Parameters
    ----------
    order : int
        Delta sigma modulator order
    h_ir : array_like of reals
        filter impulse response
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
        provide extended output, default is True
    cvxpy_xxx : various type, optional
        Parameters prefixed by ``cvxpy_`` are passed to the ``cvxpy``
        optimizer. Allowed options are:

        ``cvxpy_maxiters``
            Maximum number of iterations (defaults to 100)
        ``cvxpy_abstol``
            Absolute accuracy (defaults to 1e-7)
        ``cvxpy_reltol``
            Relative accuracy (defaults to 1e-6)
        ``cvxpy_feastol``
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break ``cvxpy`` in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.
    quad_xxx : various type
        Parameters prefixed by ``quad_`` are passed to the ``quad``
        function that is used internally as an integrator. Allowed options
        are ``quad_epsabs``, ``quad_epsrel``, ``quad_limit``, ``quad_points``.
        Do not use other options since they could break the integrator in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    See Also
    --------
    scipy.integrate.quad : integrator used internally.
        For the meaning of the integrator parameters.

    Check also the documentation of ``cvxopt`` for further information.

    Notes
    -----
    Since this function internally uses ``synthesize_ntf_from_filter``,
    the latter default parameters may also affect its behavior.
    """
    warn('Deprecated function synthesize_ntf_from_filter_imp.\n'
         'Will be removed shortly.\n'
         'Use synthesize_ntf_from_filter instead.', PyDsmWarning)
    # Manage optional parameters
    opts = synthesize_ntf_from_filter_imp.default_options.copy()
    opts.update(options)
    # Do the computation
    return synthesize_ntf_from_filter(order, h_ir, 'imp', H_inf, normalize,
                                      **opts)

synthesize_ntf_from_filter_imp.default_options = {'quad_epsabs': 1E-14,
                                                  'quad_epsrel': 1E-9}


def synthesize_ntf_from_filter_mag(order, h_mag, H_inf=1.5, normalize="auto",
                                   **options):
    u"""Synthesize a FIR NTF based on the ΔΣM output filter magnitude response.

    The ΔΣ modulator NTF is designed after the magnitude response of the
    filter in charge of removing the quantization noise

    Parameters
    ----------
    order : int
        Delta sigma modulator order
    h_mag : callable
        filter magnitude response (argument normalized in [0, 0.5])
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
        provide extended output, default is True
    cvxpy_xxx : various type, optional
        Parameters prefixed by ``cvxpy_`` are passed to the ``cvxpy``
        optimizer. Allowed options are:

        ``cvxpy_maxiters``
            Maximum number of iterations (defaults to 100)
        ``cvxpy_abstol``
            Absolute accuracy (defaults to 1e-7)
        ``cvxpy_reltol``
            Relative accuracy (defaults to 1e-6)
        ``cvxpy_feastol``
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break ``cvxpy`` in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.
    quad_xxx : various type
        Parameters prefixed by ``quad_`` are passed to the ``quad``
        function that is used internally as an integrator. Allowed options
        are ``quad_epsabs``, ``quad_epsrel``, ``quad_limit``, ``quad_points``.
        Do not use other options since they could break the integrator in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    See Also
    --------
    scipy.integrate.quad : integrator used internally.
        For the meaning of the integrator parameters.

    Check also the documentation of ``cvxopt`` for further information.

    Notes
    -----
    Since this function internally uses ``synthesize_ntf_from_filter``,
    the latter default parameters may also affect its behavior.
    """
    warn('Deprecated function synthesize_ntf_from_filter_mag.\n'
         'Will be removed shortly.\n'
         'Use synthesize_ntf_from_filter instead.', PyDsmWarning)
    # Manage optional parameters
    opts = synthesize_ntf_from_filter_mag.default_options.copy()
    opts.update(options)
    # Do the computation
    return synthesize_ntf_from_filter(order, h_mag, 'mag', H_inf, normalize,
                                      **opts)

synthesize_ntf_from_filter_mag.default_options = {'quad_epsabs': 1E-14,
                                                  'quad_epsrel': 1E-9}


def synthesize_ntf_from_filter(order, F, F_type='zpk', H_inf=1.5,
                               normalize="auto", **options):
    u"""Synthesize a FIR NTF based on the ΔΣM output filter.

    The ΔΣ modulator NTF is designed after a specification of the
    filter in charge of removing the quantization noise

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
        provide extended output, default is True
    cvxpy_xxx : various type, optional
        Parameters prefixed by ``cvxpy_`` are passed to the ``cvxpy``
        optimizer. Allowed options are:

        ``cvxpy_maxiters``
            Maximum number of iterations (defaults to 100)
        ``cvxpy_abstol``
            Absolute accuracy (defaults to 1e-7)
        ``cvxpy_reltol``
            Relative accuracy (defaults to 1e-6)
        ``cvxpy_feastol``
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break ``cvxpy`` in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.
    quad_xxx : various type
        Parameters prefixed by ``quad_`` are passed to the ``quad``
        function that is used internally as an integrator. Allowed options
        are ``quad_epsabs``, ``quad_epsrel``, ``quad_limit``, ``quad_points``.
        Do not use other options since they could break the integrator in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    See Also
    --------
    scipy.integrate.quad : integrator used internally.
        For the meaning of the integrator parameters.

    Check also the documentation of ``cvxopt`` for further information.

    Notes
    -----
    Since this function internally uses ``q0_from_filter``, the latter
    default parameters may also affect its behavior.

    Since this function internally uses ``synthesize_ntf_from_q0``, the latter
    default parameters may also affect its behavior.
    """
    # Manage optional parameters
    opts = synthesize_ntf_from_filter.default_options.copy()
    opts.update(options)
    # Do the computation
    q0 = q0_from_filter(order, F, F_type, **opts)
    return synthesize_ntf_from_q0(q0, H_inf, normalize, **opts)

synthesize_ntf_from_filter.default_options = {'quad_epsabs': 1E-14,
                                              'quad_epsrel': 1E-9}
