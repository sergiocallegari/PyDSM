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

import numpy as np
import scipy as sp
__import__("scipy.linalg")
from ._fir_weighting import q0_from_noise_weighting
import cvxpy_tinoco

__all__ = ["synthesize_ntf_from_noise_weighting",
           "synthesize_ntf_from_q0"]


def synthesize_ntf_from_q0(q0, H_inf=1.5, normalize="auto", **options):
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
    opts = synthesize_ntf_from_q0.default_options.copy()
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
    Q = cvxpy_tinoco.matrix(sp.linalg.toeplitz(q0[0:-1]))
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

synthesize_ntf_from_q0.default_options = {'cvxpy_maxiters': 100,
                                          'cvxpy_abstol': 1e-7,
                                          'cvxpy_reltol': 1e-6,
                                          'cvxpy_feastol': 1e-6,
                                          'show_progress': True}


def synthesize_ntf_from_noise_weighting(order, noise_weighting, H_inf=1.5,
                                        normalize="auto", **options):
    u"""Synthesize a FIR NTF based on a noise weighting function.

    The ΔΣ modulator NTF is designed after a noise weigthing function stating
    how expensive noise is at the various frequencies.

    Parameters
    ----------
    order : int
        Delta sigma modulator order
    noise_weighting : callable
        noise weighting function (argument normalized in [0, 0.5])
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
    """
    # Manage optional parameters
    opts = synthesize_ntf_from_noise_weighting.default_options.copy()
    opts.update(options)
    # Do the computation
    q0 = q0_from_noise_weighting(order, noise_weighting, **opts)
    return synthesize_ntf_from_q0(q0, H_inf, normalize, **opts)

synthesize_ntf_from_noise_weighting.default_options =\
    q0_from_noise_weighting.default_options.copy()
synthesize_ntf_from_noise_weighting.default_options.update(
    synthesize_ntf_from_q0.default_options)
