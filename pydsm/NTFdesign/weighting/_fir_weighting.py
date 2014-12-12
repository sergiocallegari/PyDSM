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
import scipy.linalg as la
from ...ft import idtft_hermitian
from ...delsig import evalTF, padr
import cvxpy_tinoco
from warnings import warn
from ...exceptions import PyDsmDeprecationWarning
from ...utilities import split_options, strip_options

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
            h = wi
            wn[i] = lambda f: np.abs(evalTF(h, np.exp(2j*np.pi*f)))**2
        else:
            wn[i] = wi
    return lambda f: np.prod([w(f) for w in wn], axis=0)


def q0_weighting(P, w, **options):
    """Compute Q matrix from a noise weighting function or filter

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
    quad_xxx : various type
        Parameters prefixed by ``quad_`` are passed to the ``quad``
        function that is used internally as an integrator. Allowed options
        are ``quad_epsabs``, ``quad_epsrel``, ``quad_limit``, ``quad_points``.
        Do not use other options since they could break the integrator in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    See Also
    --------
    scipy.integrate.quad : For the meaning of the integrator parameters.

    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.

    """
    # Manage parameters
    if type(w) is tuple and 2 <= len(w) <= 3:
        h = w
        w = lambda f: np.abs(evalTF(h, np.exp(2j*np.pi*f)))**2
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

    Notes
    -----
    Check the documentation of ``cvxopt`` for further information.
    """
    # Manage optional parameters
    opts = ntf_hybrid_from_q0.default_options.copy()
    opts.update(options)
    o = split_options(opts, ['cvxpy_'], ['show_progress'])
    quiet = not o.get('show_progress', True)
    # Do the computation
    order = q0.shape[0]-1
    poles = np.asarray(poles).reshape(-1)
    if poles.shape[0] > order:
        raise ValueError('Too many poles provided')
    poles = padr(poles, order, 0)
    # Get denominator coefficients from a_1 to a_order (a_0 is 1)
    ar = np.poly(poles)[1:].real
    if normalize == 'auto':
        q0 = q0/q0[0]
    elif normalize is not None:
        q0 = q0*normalize
    Q = cvxpy_tinoco.matrix(la.toeplitz(q0[0:-1]))
    L = cvxpy_tinoco.matrix(2*q0[1:])
    # Vector of numerator coefficients to be found from b_1 to b_order
    # (b_0 = 1)
    br = cvxpy_tinoco.variable(order, 1, name='br')
    # ar = cvxpy_tinoco.variable(order, 1, name='ar')
    target = cvxpy_tinoco.quad_form(br, Q) + L*br
    X = cvxpy_tinoco.variable(order, order, structure='symmetric', name='X')
    A = cvxpy_tinoco.matrix(np.eye(order, order, 1.))
    A[order-1] = -ar[::-1]
    B = cvxpy_tinoco.vstack((cvxpy_tinoco.zeros((order-1, 1)), 1.))
    C = (cvxpy_tinoco.matrix(np.eye(order, order)[:, ::-1])*br).T
    C = C-cvxpy_tinoco.matrix(ar[::-1])
    D = cvxpy_tinoco.matrix(1.)
    M1 = A.T*X
    M2 = M1*B
    M = cvxpy_tinoco.vstack((
        cvxpy_tinoco.hstack((M1*A-X, M2, C.T)),
        cvxpy_tinoco.hstack((M2.T, B.T*X*B-H_inf**2, D)),
        cvxpy_tinoco.hstack((C, D, cvxpy_tinoco.matrix(-1.)))
        ))
    constraint1 = cvxpy_tinoco.belongs(-M, cvxpy_tinoco.semidefinite_cone)
    constraint2 = cvxpy_tinoco.belongs(X, cvxpy_tinoco.semidefinite_cone)
    p = cvxpy_tinoco.program(cvxpy_tinoco.minimize(target),
                             [constraint1, constraint2])
    p.options.update(strip_options(o, 'cvxpy_'))
    p.solve(quiet)
    ntf_ir = np.hstack((1, np.asarray(br.value.T)[0]))
    return (np.roots(ntf_ir), poles, 1.)

ntf_hybrid_from_q0.default_options = {'cvxpy_maxiters': 100,
                                      'cvxpy_abstol': 1e-7,
                                      'cvxpy_reltol': 1e-6,
                                      'cvxpy_feastol': 1e-6,
                                      'show_progress': True}


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

    See Also
    --------

    Check the documentation of ``cvxopt`` for further information.
    """
    # Manage optional parameters
    opts = ntf_fir_from_q0.default_options.copy()
    opts.update(options)
    o = split_options(opts, ['cvxpy_'], ['show_progress'])
    quiet = not o.get('show_progress', True)
    # Do the computation
    order = q0.shape[0]-1
    if normalize == 'auto':
        q0 = q0/q0[0]
    elif normalize is not None:
        q0 = q0*normalize
    Q = cvxpy_tinoco.matrix(la.toeplitz(q0[0:-1]))
    L = cvxpy_tinoco.matrix(2*q0[1:])
    br = cvxpy_tinoco.variable(order, 1, name='br')
    target = cvxpy_tinoco.quad_form(br, Q) + L*br
    X = cvxpy_tinoco.variable(order, order, structure='symmetric', name='X')
    A = cvxpy_tinoco.matrix(np.eye(order, order, 1))
    B = cvxpy_tinoco.vstack((cvxpy_tinoco.zeros((order-1, 1)), 1.))
    C = (cvxpy_tinoco.matrix(np.eye(order, order)[:, ::-1])*br).T
    D = cvxpy_tinoco.matrix(1.)
    M1 = A.T*X
    M2 = M1*B
    M = cvxpy_tinoco.vstack((
        cvxpy_tinoco.hstack((M1*A-X, M2, C.T)),
        cvxpy_tinoco.hstack((M2.T, B.T*X*B-H_inf**2, D)),
        cvxpy_tinoco.hstack((C, D, cvxpy_tinoco.matrix(-1.)))
        ))
    constraint1 = cvxpy_tinoco.belongs(-M, cvxpy_tinoco.semidefinite_cone)
    constraint2 = cvxpy_tinoco.belongs(X, cvxpy_tinoco.semidefinite_cone)
    p = cvxpy_tinoco.program(cvxpy_tinoco.minimize(target),
                             [constraint1, constraint2])
    p.options.update(strip_options(o, 'cvxpy_'))
    p.solve(quiet)
    ntf_ir = np.hstack((1, np.asarray(br.value.T)[0]))
    return (np.roots(ntf_ir), np.zeros(order), 1.)

ntf_fir_from_q0.default_options = {'cvxpy_maxiters': 100,
                                   'cvxpy_abstol': 1e-7,
                                   'cvxpy_reltol': 1e-6,
                                   'cvxpy_feastol': 1e-6,
                                   'show_progress': True}


def ntf_hybrid_weighting(order, w, H_inf=1.5, poles=[],
                         normalize="auto", **options):
    u"""
    Synthesize am NTF based on a noise weighting function or filter and poles.

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
            For the meaning of the integrator parametersa.

    Notes
    -----
    Check also the documentation of ``cvxopt`` for further information.
    """
    # Manage optional parameters
    opts = ntf_fir_weighting.default_options.copy()
    opts.update(options)
    o = split_options(opts, ['quad_', 'cvxpy_'], ['show_progress'])
    # Do the computation
    poles = np.asarray(poles).reshape(-1)
    if len(poles) > 0:
        wn = mult_weightings(w, ([], poles, 1))
    else:
        wn = w
    q0 = q0_weighting(order, wn, **o['quad_'])
    return ntf_hybrid_from_q0(q0, H_inf, poles, normalize,
                              show_progress=o['show_progress'], **o['cvxpy_'])

ntf_hybrid_weighting.default_options = q0_weighting.default_options.copy()
ntf_hybrid_weighting.default_options.update(ntf_hybrid_from_q0.default_options)


def ntf_fir_weighting(order, w, H_inf=1.5,
                      normalize="auto", **options):
    u"""Synthesize a FIR NTF based on a noise weighting function or filter.

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
            For the meaning of the integrator parametersa.

    Notes
    -----
    Check also the documentation of ``cvxopt`` for further information.
    """
    # Manage optional parameters
    opts = ntf_fir_weighting.default_options.copy()
    opts.update(options)
    o = split_options(opts, ['quad_', 'cvxpy_'], ['show_progress'])
    # Do the computation
    q0 = q0_weighting(order, w, **o['quad_'])
    return ntf_fir_from_q0(q0, H_inf, normalize,
                           show_progress=o['show_progress'], **o['cvxpy_'])

ntf_fir_weighting.default_options = q0_weighting.default_options.copy()
ntf_fir_weighting.default_options.update(ntf_fir_from_q0.default_options)


# Following part is deprecated

def q0_from_noise_weighting(P, w, **options):
    warn("Function superseded by q0_weighting in "
         "NTFdesign.weighting module", PyDsmDeprecationWarning)
    return q0_weighting(P, w, **options)

q0_from_noise_weighting.__doc__ = q0_weighting.__doc__ + """
    .. deprecated:: 0.11.0
        Function has been moved to the ``NTFdesign.weighting`` module with
        name ``q0_weighting``.
    """

q0_from_noise_weighting.default_options = q0_weighting.default_options


def synthesize_ntf_from_q0(q0, H_inf=1.5, normalize="auto", **options):
    warn("Function superseded by ntf_fir_from_q0 in "
         "NTFdesign.weighting module", PyDsmDeprecationWarning)
    return ntf_fir_from_q0(q0, H_inf, normalize, **options)

synthesize_ntf_from_q0.__doc__ = ntf_fir_from_q0.__doc__ + """
    .. deprecated:: 0.11.0
        Function has been moved to the ``NTFdesign.weighting`` module with
        name ``ntf_fir_from_q0``.
    """

synthesize_ntf_from_q0.default_options = ntf_fir_from_q0.default_options


def synthesize_ntf_from_noise_weighting(order, noise_weighting, H_inf=1.5,
                                        normalize="auto", **options):
    warn("Function superseded by ntf_fir_weighting in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_fir_weighting(order, noise_weighting, H_inf,
                             normalize, **options)

synthesize_ntf_from_noise_weighting.__doc__ = ntf_fir_weighting.__doc__ + """
    .. deprecated:: 0.11.0
        Function has been moved to the ``NTFdesign`` module with
        name ``ntf_fir_weighting``.
    """

synthesize_ntf_from_noise_weighting.default_options = \
    ntf_fir_weighting.default_options
