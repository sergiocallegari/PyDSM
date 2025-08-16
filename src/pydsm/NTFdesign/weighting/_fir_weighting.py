# -*- coding: utf-8 -*-

# Copyright (c) 2013–2024, Sergio Callegari
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
from ...delsig import evalTF, padr
from warnings import warn
from ...exceptions import PyDsmDeprecationWarning
from ...utilities import digested_options
import scipy.linalg as la

__all__ = ["q0_from_noise_weighting", "q0_weighting",
           "ntf_fir_from_q0", "synthesize_ntf_from_q0",
           "ntf_fir_weighting", "synthesize_ntf_from_noise_weighting",
           "mult_weightings", "ntf_hybrid_weighting"]


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
        integrator in unexpected ways.

    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.

    Default values for the options not directly documented in the function
    call signature can be checked and updated by changing the function
    ``default_options`` attribute.

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
    return idtft_hermitian(w, np.arange(P+1), **opts)

q0_weighting.default_options = {"quad_opts": {"epsabs": 1E-14,
                                              "epsrel": 1E-9,
                                              "limit": 100,
                                              "points": None}}


def ntf_fir_from_q0(q0, H_inf=1.5, normalize="auto", **options):
    """Synthesize FIR NTF from quadratic form expressing noise weighting.

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
        provide extended output.
    fix_pos : bool, optional
        fix quadratic form for positive definiteness. Numerical noise
        may make it not positive definite leading to errors.
    modeler : string, optional
        modeling backend for the optimization problem. Currently, the
        ``cvxpy_old``, ``cvxpy`` and ``picos`` backends are supported.
        Default is ``cvxpy_old``.
    cvxpy_opts : dictionary, optional
       A dictionary of options to use with the ``cvxpy`` modeling library.
       Allowed options include:

       ``override_kktsolver`` (bool)
           Whether to override the default ``cvxopt`` kkt solver using the
           ``chol`` kkt solver.
           Leave this at the default True setting, to avoid paying a
           performance price.
       ``solver`` (string)
           The solver backend to use. Either `cvxopt` or `scs`

    cvxopt_opts : dict, optional
        A dictionary of options for the ``cvxopt`` optimizer.
        Allowed options include:

        ``maxiters`` (int)
            Maximum number of iterations
        ``abstol`` (real)
            Absolute accuracy
        ``reltol`` (real)
            Relative accuracy
        ``feastol`` (real)
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break ``CVXOPT`` in
        unexpected ways. These options can be passed when using the
        ``cvxpy_old`` modeler, the ``picos`` modeler or the ``cvxpy`` modeler
        with the ``cvxopt`` backend.

    scs_opts : dict, optional
        A dictionary of options for the ``scs`` optimizer.  Allowed options
        include:

        ``max_iters`` (int)
            Maximum number of iterations
        ``eps`` (real)
            Convergence tolerance
        ``alpha`` (real)
            Relaxation parameter
        ``normalize`` (bool)
            Whether to precondition data matrices
        ``use_indirect`` (bool)
            Whether to use indirect solver for KKT sytem (instead of direct)

       Do not use other options since they could break ``scs`` in
       unexpected ways. These options can be passed when using the
       ``cvxpy`` modeler with the ``scs`` backend.

    Notes
    -----
    Default values for the options not directly documented in the
    function call signature can be checked and updated by changing the
    function ``default_options`` attribute. For more information on
    the specific options taken by the ``CVXOPT`` and ``SCS``
    optimizer, look at the corresponding documentation. Similarly, for
    more information about the ``cvxpy`` modeler and optimization
    frontend, look at the ``cvxpy`` documentation.
    """
    # Manage optional parameters
    opts = digested_options(
        options, ntf_fir_from_q0.default_options,
        ['show_progress', 'fix_pos', 'modeler'], [], False)
    dig_opts = {'show_progress': opts['show_progress'],
                'cvxpy_opts': {},
                'cvxpy_tdr_opts': {},
                'picos_opts': {}}
    if opts['modeler'] == 'cvxpy':
        opts.update(digested_options(
            options, ntf_fir_from_q0.default_options,
            [], ['cvxpy_opts'], False))
        if opts['cvxpy_opts']['solver'] == 'cvxopt':
            dig_opts['cvxpy_opts'].update(digested_options(
                options, ntf_fir_from_q0.default_options,
                [], ['cvxopt_opts'], False)['cvxopt_opts'])
            if opts['cvxpy_opts'].get('override_kktsolver', True):
                dig_opts['cvxpy_opts']['kktsolver'] = 'chol'
        elif opts['cvxpy_opts']['solver'] == 'scs':
            dig_opts['cvxpy_opts'].update(digested_options(
                options, ntf_fir_from_q0.default_options,
                [], ['scs_opts'], False)['scs_opts'])
        opts['cvxpy_opts'].pop('override_kktsolver')
        dig_opts['cvxpy_opts'].update(opts['cvxpy_opts'])
        from ._fir_weighting_cvxpy import (
            ntf_fir_from_digested as _ntf_fir_from_digested)
    elif opts['modeler'] == 'cvxpy_old':
        dig_opts['cvxpy_tdr_opts'].update(digested_options(
            options, ntf_fir_from_q0.default_options,
            [], ['cvxopt_opts'], False)['cvxopt_opts'])
        from ._fir_weighting_cvxpy_tdr import (
            ntf_fir_from_digested as _ntf_fir_from_digested)
    elif opts['modeler'] == 'picos':
        dig_opts['picos_opts'].update(digested_options(
            options, ntf_fir_from_q0.default_options,
            [], ['cvxopt_opts'], False)['cvxopt_opts'])
        from ._fir_weighting_picos import (
            ntf_fir_from_digested as _ntf_fir_from_digested)
    else:
        raise ValueError('Unsupported modeling backend {}'.format(
            opts['modeler']))
    digested_options(options, {})
    # Do the computation
    if normalize == 'auto':
        q0 = q0/q0[0]
    elif normalize is not None:
        q0 = q0*normalize
    order = q0.shape[0]-1
    Q = la.toeplitz(q0)
    d, v = np.linalg.eigh(Q)
    if opts['fix_pos']:
        d = d/np.max(d)
        d[d < 0] = 0.
    Qs = v.dot(np.diag(np.sqrt(d))).dot(v.T)
    A = np.eye(order, order, 1)
    C = np.zeros((1, order))
    ntf_ir = _ntf_fir_from_digested(Qs, A, C, H_inf, **dig_opts)
    return (np.roots(ntf_ir), np.zeros(order), 1.)


ntf_fir_from_q0.default_options = {"modeler": "cvxpy_old",
                                   "cvxpy_opts": {'override_kktsolver': True,
                                                  'solver': 'cvxopt'},
                                   "cvxopt_opts": {'maxiters': 100,
                                                   'abstol': 1e-7,
                                                   'reltol': 1e-6,
                                                   'feastol': 1e-6},
                                   'scs_opts': {'max_iters': 250000,
                                                'eps': 1e-10,
                                                'alpha': 1.8,
                                                'normalize': True,
                                                'use_indirect': False,
                                                'acceleration_lookback': 0},
                                   'show_progress': True,
                                   'fix_pos': True}


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
        provide extended output.
    fix_pos : bool, optional
        fix quadratic form for positive definiteness. Numerical noise
        may make it not positive definite leading to errors.
    modeler : string, optional
        modeling backend for the optimization problem. Currently, the
        ``cvxpy_old``, ``cvxpy`` and ``picos`` backends are supported.
        Default is ``cvxpy_old``.
    cvxpy_opts : dictionary, optional
       A dictionary of options to use with the ``cvxpy`` modeling library.
       Allowed options include:

       ``override_kktsolver`` (bool)
           Whether to override the default ``cvxopt`` kkt solver using the
           ``chol`` kkt solver.
           Leave this at the default True setting, to avoid paying a
           performance price.
       ``solver`` (string)
           The solver backend to use. Either `cvxopt` or `scs`

    cvxopt_opts : dict, optional
        A dictionary of options for the ``cvxopt`` optimizer.
        Allowed options include:

        ``maxiters`` (int)
            Maximum number of iterations
        ``abstol`` (real)
            Absolute accuracy
        ``reltol`` (real)
            Relative accuracy
        ``feastol`` (real)
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break ``cvxopt`` in
        unexpected ways. These options can be passed when using the
        ``cvxpy_old`` modeler, the ``picos`` modeler or the ``cvxpy`` modeler
        with the ``cvxopt`` backend.
    scs_opts : dict, optional
        A dictionary of options for the ``scs`` optimizer.  Allowed options
        include:

        ``max_iters`` (int)
            Maximum number of iterations
        ``eps`` (real)
            Convergence tolerance
        ``alpha`` (real)
            Relaxation parameter
        ``normalize`` (bool)
            Whether to precondition data matrices
        ``use_indirect`` (bool)
            Whether to use indirect solver for KKT sytem (instead of direct)

       Do not use other options since they could break ``scs`` in
       unexpected ways. These options can be passed when using the
       ``cvxpy`` modeler with the ``scs`` backend.
    quad_opts : dictionary, optional
        Parameters to be passed to the ``quad`` function used internally as
        an integrator. Allowed options are ``epsabs``, ``epsrel``, ``limit``,
        ``points``. Do not use other options since they could break the
        integrator in unexpected ways.

    Notes
    -----
    Default values for the options not directly documented in the function
    call signature can be checked and updated by changing the function
    ``default_options`` attribute.

    The internal operation of this function is described in [1]_.

    .. [1] Sergio Callegari, Federico Bizzarri *“Noise Weighting in the
       Design of ΔΣ Modulators (with a Psychoacoustic Coder as an
       Example),”* IEEE Transactions on Circuits and Systems - Part II:
       Express Briefs, Vol. 60, N. 11, pp. 756-760. Nov. 2013. DOI:
       `10.1109/TCSII.2013.2281892
       <http://dx.doi.org/10.1109/TCSII.2013.2281892>`__. Pre-print
       available on `arXiv <http://arxiv.org/abs/1309.6151>`__.

    For the specific options taken by the modeler and optimization
    frontend ``cvxpy`` as well as the optimizers ``CVXOPT`` and
    ``scs``, look at the documentation of the individual packages.

    See Also
    --------
    scipy.integrate.quad : for the meaning of the integrator parameters
    """
    # Manage optional parameters
    opts1 = digested_options(options, ntf_fir_weighting.default_options,
                             [], ['quad_opts'], False)
    opts2 = digested_options(
        options, ntf_fir_weighting.default_options,
        ['show_progress', 'fix_pos', 'modeler'], [], False)
    if opts2['modeler'] == 'cvxpy':
        opts2.update(digested_options(
            options, ntf_fir_weighting.default_options,
            [], ['cvxpy_opts'], False))
        if opts2['cvxpy_opts']['solver'] == 'cvxopt':
            opts2.update(digested_options(
                options, ntf_fir_weighting.default_options,
                [], ['cvxopt_opts'], False))
        elif opts2['cvxpy_opts']['solver'] == 'scs':
            opts2.update(digested_options(
                options, ntf_fir_weighting.default_options,
                [], ['scs_opts'], False))
    elif opts2['modeler'] == 'cvxpy_old' or opts2['modeler'] == 'picos':
            opts2.update(digested_options(
                options, ntf_fir_weighting.default_options,
                [], ['cvxopt_opts'], False))
    else:
        raise ValueError('Unsupported modeling backend {}'.format(
            opts2['modeler']))
    digested_options(options, {})
    # Do the computation
    q0 = q0_weighting(order, w, **opts1)
    return ntf_fir_from_q0(q0, H_inf, normalize, **opts2)

ntf_fir_weighting.default_options = q0_weighting.default_options.copy()
ntf_fir_weighting.default_options.update(ntf_fir_from_q0.default_options)


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
        provide extended output.
    fix_pos : bool, optional
        fix quadratic form for positive definiteness. Numerical noise
        may make it not positive definite leading to errors.
    modeler : string, optional
        modeling backend for the optimization problem. Currently, the
        ``cvxpy_old``, ``cvxpy`` and ``picos`` backends are supported.
        Default is ``cvxpy_old``.
    cvxpy_opts : dictionary, optional
       A dictionary of options to use with the ``cvxpy`` modeling library.
       Allowed options include:

       ``override_kktsolver`` (bool)
           Whether to override the default ``cvxopt`` kkt solver using the
           ``chol`` kkt solver.
           Leave this at the default True setting, to avoid paying a
           performance price.
       ``solver`` (string)
           The solver backend to use. Either `cvxopt` or `scs`

    cvxopt_opts : dict, optional
        A dictionary of options for the ``cvxopt`` optimizer.
        Allowed options include:

        ``maxiters`` (int)
            Maximum number of iterations
        ``abstol`` (real)
            Absolute accuracy
        ``reltol`` (real)
            Relative accuracy
        ``feastol`` (real)
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break ``cvxopt`` in
        unexpected ways. These options can be passed when using the
        ``cvxpy_old`` modeler, the ``picos`` modeler or the ``cvxpy`` modeler
        with the ``cvxopt`` backend.
    scs_opts : dict, optional
        A dictionary of options for the ``scs`` optimizer.  Allowed options
        include:

        ``max_iters`` (int)
            Maximum number of iterations
        ``eps`` (real)
            Convergence tolerance
        ``alpha`` (real)
            Relaxation parameter
        ``normalize`` (bool)
            Whether to precondition data matrices
        ``use_indirect`` (bool)
            Whether to use indirect solver for KKT sytem (instead of direct)

       Do not use other options since they could break ``scs`` in
       unexpected ways. These options can be passed when using the
       ``cvxpy`` modeler with the ``scs`` backend.
    quad_opts : dictionary, optional
        Parameters to be passed to the ``quad`` function used internally as
        an integrator. Allowed options are ``epsabs``, ``epsrel``, ``limit``,
        ``points``. Do not use other options since they could break the
        integrator in unexpected ways.

    Notes
    -----
    Default values for the options not directly documented in the function
    call signature can be checked and updated by changing the function
    ``default_options`` attribute.

    The internal operation of this function is described in [1]_.

    .. [1] Sergio Callegari, Federico Bizzarri "Optimal Design of the Noise
       Transfer Function of ΔΣ Modulators: IIR strategies, FIR strategies,
       FIR strategies with Preassigned Poles," Signal Processing, Elsevier,
       2015. DOI:  10.1016/j.sigpro.2015.02.001.

    For the specific options taken by the modeler and optimization
    frontend ``cvxpy`` as well as the optimizers ``CVXOPT`` and
    ``scs``, look at the documentation of the individual packages.

    See Also
    --------
    scipy.integrate.quad : for the meaning of the integrator parameters
    """
    # Manage optional parameters
    opts1 = digested_options(options, ntf_hybrid_weighting.default_options,
                             [], ['quad_opts'], False)
    opts2 = digested_options(
        options, ntf_hybrid_weighting.default_options,
        ['show_progress', 'fix_pos', 'modeler'], [], False)
    dig_opts = {'show_progress': opts2['show_progress'],
                'cvxpy_opts': {},
                'cvxpy_tdr_opts': {},
                'picos_opts': {}}
    if opts2['modeler'] == 'cvxpy':
        opts2.update(digested_options(
            options, ntf_fir_from_q0.default_options,
            [], ['cvxpy_opts'], False))
        if opts2['cvxpy_opts']['solver'] == 'cvxopt':
            dig_opts['cvxpy_opts'].update(digested_options(
                options, ntf_fir_from_q0.default_options,
                [], ['cvxopt_opts'], False)['cvxopt_opts'])
            if opts2['cvxpy_opts'].get('override_kktsolver', True):
                dig_opts['cvxpy_opts']['kktsolver'] = 'chol'
        elif opts2['cvxpy_opts']['solver'] == 'scs':
            dig_opts['cvxpy_opts'].update(digested_options(
                options, ntf_fir_from_q0.default_options,
                [], ['scs_opts'], False)['scs_opts'])
        opts2['cvxpy_opts'].pop('override_kktsolver')
        dig_opts['cvxpy_opts'].update(opts2['cvxpy_opts'])
        from ._fir_weighting_cvxpy import (
            ntf_fir_from_digested as _ntf_fir_from_digested)
    elif opts2['modeler'] == 'cvxpy_old':
        dig_opts['cvxpy_tdr_opts'].update(digested_options(
            options, ntf_fir_from_q0.default_options,
            [], ['cvxopt_opts'], False)['cvxopt_opts'])
        from ._fir_weighting_cvxpy_tdr import (
            ntf_fir_from_digested as _ntf_fir_from_digested)
    elif opts2['modeler'] == 'picos':
        dig_opts['picos_opts'].update(digested_options(
            options, ntf_fir_from_q0.default_options,
            [], ['cvxopt_opts'], False)['cvxopt_opts'])
        from ._fir_weighting_picos import (
            ntf_fir_from_digested as _ntf_fir_from_digested)
    else:
        raise ValueError('Unsupported modeling backend {}'.format(
            opts2['modeler']))
    digested_options(options, {})
    # Do the computation
    poles = np.asarray(poles).reshape(-1)
    if len(poles) > 0:
        wn = mult_weightings(w, ([], poles, 1))
    else:
        wn = w
    q0 = q0_weighting(order, wn, **opts1)
    if normalize == 'auto':
        q0 = q0/q0[0]
    elif normalize is not None:
        q0 = q0*normalize
    if poles.shape[0] > order:
        raise ValueError('Too many poles provided')
    poles = padr(poles, order, 0)
    # Get denominator coefficients from a_1 to a_order (a_0 is 1)
    ar = np.poly(poles)[1:].real
    Q = la.toeplitz(q0)
    d, v = np.linalg.eigh(Q)
    if opts2['fix_pos']:
        d = d/np.max(d)
        d[d < 0] = 0.
    Qs = v.dot(np.diag(np.sqrt(d))).dot(v.T)
    A = np.eye(order, order, 1)
    A[order-1] = -ar[::-1]
    C = -ar[::-1].reshape((1, order))
    ntf_ir = _ntf_fir_from_digested(Qs, A, C, H_inf, **dig_opts)
    return (np.roots(ntf_ir), poles, 1.)

ntf_hybrid_weighting.default_options = {"modeler": "cvxpy_old",
                                        "cvxpy_opts": {'override_kktsolver':
                                                       True,
                                                       'solver': 'cvxopt'},
                                        "cvxopt_opts": {'maxiters': 100,
                                                        'abstol': 1e-7,
                                                        'reltol': 1e-6,
                                                        'feastol': 1e-6},
                                        'scs_opts': {'max_iters': 2500,
                                                     'eps': 1e-3,
                                                     'alpha': 1.8,
                                                     'normalize': True,
                                                     'use_indirect': False},
                                        'show_progress': True,
                                        'fix_pos': True}
ntf_hybrid_weighting.default_options.update(q0_weighting.default_options)


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
