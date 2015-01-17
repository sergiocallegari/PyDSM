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
Design optimal NTF with NTF zeros as parameters

Routines to design an optimal NTF from a weighting function using
the NTF zeros as free parameters.

The optimization is practiced in two steps:

* First a matrix (symmetric, Toeplitz, Positive semidefinite) is
  extracted from the weighting function to define a quadratic cost
  function for the optimization problem. The Toeplitz-symmetric property of
  the matrix let it be fully defined from its first row only.
* Then the NTF zeros are optimized based on this cost function.
"""

from __future__ import division, print_function

import numpy as np
from ...ft import idtft_hermitian
from ...delsig import evalTF
from warnings import warn
from ...exceptions import PyDsmDeprecationWarning
from ...utilities import digested_options

__all__ = ["q0_from_noise_weighting", "q0_weighting",
           "ntf_fir_from_q0", "synthesize_ntf_from_q0",
           "ntf_fir_weighting", "synthesize_ntf_from_noise_weighting",
           "mult_weightings", "ntf_hybrid_from_q0", "ntf_hybrid_weighting"]


def mult_weightings(*ww):
    """
    Product of weighting functions.

    Returns a function that is the product of many weighting functions.

    Parameters
    ----------
    w1, w2, ... : functions or tuples
        Weighting functions.
        If an entry is a function, it is used as is.
        If it is a tuple, it interpreted as filters in ba or zpk form from
        which a weighting function is implicitly obtained.

    Returns
    -------
    w : function
        Overall weighting function
    """
    wn = [0] * len(ww)
    for i, wi in enumerate(ww):
        if type(wi) is tuple and 2 <= len(wi) <= 3:
            wn[i] = (lambda wi:
                     lambda f: np.abs(evalTF(wi, np.exp(2j*np.pi*f)))**2)(wi)
        else:
            wn[i] = wi
    return lambda f: np.prod([w(f) for w in wn], axis=0)


def q0_weighting(P, w, **options):
    """Compute Q matrix from a noise weighting function or a filter

    Parameters
    ----------
    P : int
        order of the FIR to be eventually synthesized
    w : callable with argument f in [0,1/2] or tuple
            * if function: noise weighting function
            * if filter definition as zpk or ba tuple: weighting is implicitly
              provided by the filter

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

    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.

    See Also
    --------
    scipy.integrate.quad : For the meaning of the integrator parameters.
    """
    # Manage parameters
    if type(w) is tuple and 2 <= len(w) <= 3:
        h = w
        w = lambda f: np.abs(evalTF(h, np.exp(2j*np.pi*f)))**2
    # Manage optional parameters
    opts = digested_options(options, q0_weighting.default_options,
                            [], ['quad_opts'])
    # Do the computation
    ac = lambda t: idtft_hermitian(w, t, **opts)
    return np.asarray(map(ac, np.arange(P+1)))

q0_weighting.default_options = {"quad_opts": {"epsabs": 1E-14,
                                              "epsrel": 1E-9,
                                              "limit": 100,
                                              "points": None}}


def ntf_hybrid_from_q0(q0, H_inf=1.5, poles=[], normalize="auto", **options):
    """
    Synthesize NTF from quadratic form expressing noise weighting and poles.

    Parameters
    ----------
    q0 : ndarray
        first row of the Toeplitz symmetric matrix defining the quadratic form
    H_inf : real, optional
        Max peak NTF gain, defaults to 1.5, used to enforce the Lee criterion
    poles : array_like
        List of pre-assigned NTF poles. Must be no longer than length of q0
        minus 1.
    normalize : string or real, optional
        Normalization to apply to the quadratic form used in the NTF
        selection. Defaults to 'auto' which means setting the top left entry
        in the matrix Q defining the quadratic form to 1.

    Returns
    -------
    ntf : ndarray
        NTF in zpk form

    Other parameters
    ----------------
    show_progress : bool, optional
        provide extended output, default is True and can be updated by
        changing the function ``default_options`` attribute.
    fix_pos : bool, optional
        fix quadratic form for positive definiteness. Numerical noise
        may make it not positive definite leading to errors. Default is True
        and can be updated by changing the function ``default_options``
        attribute.
    modeler : string
        modeling backend for the optimization problem. Currently, only the
        ``cvxpy_old`` backend is supported.
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

        Do not use other options since they could break ``cvxpy`` in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    See also
    --------
    cvxopt : for the optimizer parameters.
    """
    # Manage optional parameters
    opts = digested_options(options, ntf_hybrid_from_q0.default_options,
                            ['show_progress', 'fix_pos', 'modeler'],
                            ['cvxopt_opts'])
    if normalize == 'auto':
        q0 = q0/q0[0]
    elif normalize is not None:
        q0 = q0*normalize
    if opts['modeler'] == 'cvxpy_old':
        from ._fir_weighting_tinoco import (
            ntf_hybrid_from_q0 as _ntf_hybrid_from_q0)
        return _ntf_hybrid_from_q0(q0, H_inf, poles, **opts)
    else:
        raise ValueError("Unsupported modeling backend")


ntf_hybrid_from_q0.default_options = {"modeler": "cvxpy_old",
                                      "cvxopt_opts": {'maxiters': 100,
                                                      'abstol': 1e-7,
                                                      'reltol': 1e-6,
                                                      'feastol': 1e-6},
                                      'show_progress': True,
                                      'fix_pos': True}


def ntf_fir_from_q0(q0, H_inf=1.5, normalize="auto", **options):
    """
    Synthesize FIR NTF from quadratic form expressing noise weighting.

    Parameters
    ----------
    q0 : ndarray
        first row of the Toeplitz symmetric matrix defining the quadratic form
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
    fix_pos : bool, optional
        fix quadratic form for positive definiteness. Numerical noise
        may make it not positive definite leading to errors. Default is True
        and can be updated by changing the function ``default_options``
        attribute.
    modeler : string
        modeling backend for the optimization problem. Currently, the
        ``cvxpy_old`` and ``cvxpy`` backends are supported. Default is
        ``cvxpy_old``.
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

        Do not use other options since they could break ``cvxpy`` in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    See Also
    --------
    cvxopt : for the optimizer parameters
    """
    # Manage optional parameters
    opts = digested_options(options, ntf_fir_from_q0.default_options,
                            ['show_progress', 'fix_pos', 'modeler'],
                            ['cvxopt_opts'])
    # Do the computation
    if normalize == 'auto':
        q0 = q0/q0[0]
    elif normalize is not None:
        q0 = q0*normalize
    if opts['modeler'] == 'cvxpy_old':
        from ._fir_weighting_tinoco import ntf_fir_from_q0 as _ntf_fir_from_q0
    elif opts['modeler'] == 'cvxpy':
        from ._fir_weighting_cvxpy import ntf_fir_from_q0 as _ntf_fir_from_q0
    else:
        raise ValueError("Unsupported modeling backend")
    return _ntf_fir_from_q0(q0, H_inf, **opts)


ntf_fir_from_q0.default_options = {"modeler": "cvxpy_old",
                                   "cvxopt_opts": {'maxiters': 100,
                                                   'abstol': 1e-7,
                                                   'reltol': 1e-6,
                                                   'feastol': 1e-6},
                                   'show_progress': True,
                                   'fix_pos': True}


def ntf_hybrid_weighting(order, w, H_inf=1.5, poles=[],
                         normalize="auto", **options):
    u"""
    Synthesize NTF based on noise weighting function or filter plus poles.

    The ΔΣ modulator NTF is designed after a noise weigthing function stating
    how expensive noise is at the various frequencies.

    Parameters
    ----------
    order : int
        Delta sigma modulator order
    w : callable with argument f in [0,1/2] or tuple
            * if function: noise weighting function
            * if filter definition as zpk or ba tuple: weighting is implicitly
              provided by the filter
    H_inf : real, optional
        Max peak NTF gain, defaults to 1.5, used to enforce the Lee criterion
    poles : array_like
        List of pre-assigned NTF poles. Must be no longer than ``order``.
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
    fix_pos : bool, optional
        fix quadratic form for positive definiteness. Numerical noise
        may make it not positive definite leading to errors. Default is True
        and can be updated by changing the function ``default_options``
        attribute.
    modeler : string
        modeling backend for the optimization problem. Currently, only the
        ``cvxpy_old`` backend is supported.
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

        Do not use other options since they could break ``cvxpy`` in
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
    scipy.integrate.quad : for the meaning of the integrator parameters
    cvxopt : for the optimizer parameters
    """
    # Manage optional parameters
    opts1 = digested_options(options, ntf_hybrid_weighting.default_options,
                             [], ['quad_opts'], False)
    opts2 = digested_options(options, ntf_hybrid_weighting.default_options,
                             ['show_progress', 'fix_pos', 'modeler'],
                             ['cvxopt_opts'])
    # Do the computation
    poles = np.asarray(poles).reshape(-1)
    if len(poles) > 0:
        wn = mult_weightings(w, ([], poles, 1))
    else:
        wn = w
    q0 = q0_weighting(order, wn, **opts1)
    return ntf_hybrid_from_q0(q0, H_inf, poles, normalize, **opts2)

ntf_hybrid_weighting.default_options = q0_weighting.default_options.copy()
ntf_hybrid_weighting.default_options.update(ntf_hybrid_from_q0.default_options)


def ntf_fir_weighting(order, w, H_inf=1.5,
                      normalize="auto", **options):
    u"""Synthesize FIR NTF based on a noise weighting function or a filter.

    The ΔΣ modulator NTF is designed after a noise weigthing function stating
    how expensive noise is at the various frequencies.

    Parameters
    ----------
    order : int
        Delta sigma modulator order
    w : callable with argument f in [0,1/2] or tuple
            * if function: noise weighting function
            * if filter definition as zpk or ba tuple: weighting is implicitly
              provided by the filter
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
    fix_pos : bool, optional
        fix quadratic form for positive definiteness. Numerical noise
        may make it not positive definite leading to errors. Default is True
        and can be updated by changing the function ``default_options``
        attribute.
    modeler : string
        modeling backend for the optimization problem. Currently, the
        ``cvxpy_old`` and ``cvxpy`` backends are supported. Default is
        ``cvxpy_old``.
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

        Do not use other options since they could break ``cvxpy`` in
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
    scipy.integrate.quad : for the meaning of the integrator parameters
    cvxopt : for the optimizer parameters
    """
    # Manage optional parameters
    opts1 = digested_options(options, ntf_fir_weighting.default_options,
                             [], ['quad_opts'], False)
    opts2 = digested_options(options, ntf_fir_weighting.default_options,
                             ['show_progress', 'fix_pos', 'modeler'],
                             ['cvxopt_opts'])
    # Do the computation
    q0 = q0_weighting(order, w, **opts1)
    return ntf_fir_from_q0(q0, H_inf, normalize, **opts2)

ntf_fir_weighting.default_options = q0_weighting.default_options.copy()
ntf_fir_weighting.default_options.update(ntf_fir_from_q0.default_options)


# Following part is deprecated

def q0_from_noise_weighting(P, w, **options):
    """
    Alias of :func:`q0_weighting`

    .. deprecated:: 0.11.0
        Function has been renamed :func:`q0_weighting`.
    """
    warn("Function superseded by q0_weighting in "
         "NTFdesign.weighting module", PyDsmDeprecationWarning)
    return q0_weighting(P, w, **options)


q0_from_noise_weighting.default_options = q0_weighting.default_options


def synthesize_ntf_from_q0(q0, H_inf=1.5, normalize="auto", **options):
    """
    Alias of :func:`ntf_fir_from_q0`

    .. deprecated:: 0.11.0
        Function has been renamed :func:`ntf_fir_from_q0`.
    """
    warn("Function superseded by ntf_fir_from_q0 in "
         "NTFdesign.weighting module", PyDsmDeprecationWarning)
    return ntf_fir_from_q0(q0, H_inf, normalize, **options)


def synthesize_ntf_from_noise_weighting(order, noise_weighting, H_inf=1.5,
                                        normalize="auto", **options):
    """
    Alias of :func:`ntf_fir_weighting`

    .. deprecated:: 0.11.0
        Function has been renamed :func:`ntf_fir_weighting`.
    """
    warn("Function superseded by ntf_fir_weighting in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_fir_weighting(order, noise_weighting, H_inf,
                             normalize, **options)
