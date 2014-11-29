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
from ...delsig import evalTF
import cvxpy_tinoco
from warnings import warn
from ...exceptions import PyDsmDeprecationWarning

__all__ = ["q0_from_noise_weighting", "q0_weighting",
           "ntf_fir_from_q0", "synthesize_ntf_from_q0",
           "ntf_fir_weighting", "synthesize_ntf_from_noise_weighting"]


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

    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.

    See Also
    --------
    scipy.integrate.quad : integrator used internally.
        For the meaning of the integrator parameters.
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


def ntf_fir_from_q0(q0, H_inf=1.5, normalize="auto", **options):
    """Synthesize FIR NTF from quadratic form defining quantization noise gain.

    Parameters
    ----------
    q0 : array_like
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
    cvxpy_opts = {k[6:]: v for k, v in opts.iteritems()
                  if k.startswith('cvxpy_')}
    if 'show_progress' in opts:
        quiet = not opts['show_progress']
    else:
        quiet = False
    # Do the computation
    order = q0.shape[0]-1
    if normalize == 'auto':
        q0 = q0/q0[0]
    elif normalize is not None:
        q0 = q0*normalize
    Q = cvxpy_tinoco.matrix(la.toeplitz(q0[0:-1]))
    L = cvxpy_tinoco.matrix(2*q0[1:])
    ar = cvxpy_tinoco.variable(order, 1, name='ar')
    target = cvxpy_tinoco.quad_form(ar, Q) + L*ar
    X = cvxpy_tinoco.variable(order, order, structure='symmetric', name='X')
    A = cvxpy_tinoco.matrix(np.eye(order, order, 1))
    B = cvxpy_tinoco.matrix(np.vstack((np.zeros((order-1, 1)), 1)))
    C = (cvxpy_tinoco.matrix(np.eye(order, order)[:, ::-1])*ar).T
    D = cvxpy_tinoco.matrix([[1]])
    M1 = A.T*X
    M2 = M1*B
    M = cvxpy_tinoco.vstack((
        cvxpy_tinoco.hstack((M1*A-X, M2, C.T)),
        cvxpy_tinoco.hstack((M2.T, B.T*X*B-H_inf**2, D)),
        cvxpy_tinoco.hstack((C, D, cvxpy_tinoco.matrix([[-1]])))
        ))
    constraint1 = cvxpy_tinoco.belongs(-M, cvxpy_tinoco.semidefinite_cone)
    constraint2 = cvxpy_tinoco.belongs(X, cvxpy_tinoco.semidefinite_cone)
    p = cvxpy_tinoco.program(cvxpy_tinoco.minimize(target),
                             [constraint1, constraint2])
    p.options.update(cvxpy_opts)
    p.solve(quiet)
    ntf_ir = np.hstack((1, np.asarray(ar.value.T)[0]))
    return (np.roots(ntf_ir), np.zeros(order), 1.)

ntf_fir_from_q0.default_options = {'cvxpy_maxiters': 100,
                                   'cvxpy_abstol': 1e-7,
                                   'cvxpy_reltol': 1e-6,
                                   'cvxpy_feastol': 1e-6,
                                   'show_progress': True}


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
    # Do the computation
    q0 = q0_weighting(order, w, **opts)
    return ntf_fir_from_q0(q0, H_inf, normalize, **opts)

ntf_fir_weighting.default_options = q0_weighting.default_options.copy()
ntf_fir_weighting.default_options.update(ntf_fir_from_q0.default_options)


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


def synthesize_ntf_from_q0(q0, H_inf=1.5, normalize="auto", **options):
    warn("Function superseded by ntf_fir_from_q0 in "
         "NTFdesign.weighting module", PyDsmDeprecationWarning)
    return ntf_fir_from_q0(q0, H_inf, normalize, **options)

synthesize_ntf_from_q0.__doc__ = ntf_fir_from_q0.__doc__ + """
    .. deprecated:: 0.11.0
    Function has been moved to the ``NTFdesign.weighting`` module with name
    ``ntf_fir_from_q0``.
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
    Function has been moved to the ``NTFdesign`` module with name
    ``ntf_fir_weighting``.
    """

synthesize_ntf_from_noise_weighting.default_options = \
    ntf_fir_weighting.default_options
