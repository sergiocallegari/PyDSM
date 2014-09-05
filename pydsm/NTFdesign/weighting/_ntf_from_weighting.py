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

import numpy as np
import scipy as sp
__import__("scipy.linalg")
from ._q0_from_weighting import q0_from_noise_weighting
import cvxpy

__all__=["synthesize_ntf_from_noise_weighting",
         "synthesize_ntf_from_q0"]


def synthesize_ntf_from_noise_weighting(order, noise_weighting, H_inf=1.5,
                                        normalize="auto", options={}):
    u"""
    Synthesize a FIR NTF based on a noise weighting function.

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
    options : dict, optional
        parameters for the SDP optimizer, see the documentation of `cvxpy`.
        This includes 'show_progress' (default True).

    Returns
    -------
    ntf : ndarray
        FIR NTF in zpk form

    Notes
    -----
    The computation of the NTF from the noise weighting involves computing
    an integral on the noise weighting function. To control the integration
    parameters, do not use this function. Rather, first compute a vector
    q0 with `q0_from_noise_weighting` (which lets the integrator params be
    specified), then use `synthesize_ntf_from_q0`.
    """
    q0=q0_from_noise_weighting(order, noise_weighting)
    return synthesize_ntf_from_q0(q0, H_inf, normalize, options)


def synthesize_ntf_from_q0(q0, H_inf=1.5, normalize="auto",
                          options={}):
    """
    Synthesize a FIR NTF from quadratic form defining quantization noise gain.

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
    options : dict, optional
        parameters for the SDP optimizer, see the documentation of `cvxpy`.
        This includes 'show_progress' (default True).

    Returns
    -------
    ntf : ndarray
        FIR NTF in zpk form
    """
    order=q0.shape[0]-1
    if normalize=='auto':
        q0=q0/q0[0]
    elif normalize!=None:
        q0=q0*normalize
    Q=cvxpy.matrix(sp.linalg.toeplitz(q0[0:-1]))
    L=cvxpy.matrix(2*q0[1:])
    ar=cvxpy.variable(order,1, name='ar')
    target=cvxpy.quad_form(ar,Q) + L*ar
    X=cvxpy.variable(order, order, structure='symmetric', name='X')
    A=cvxpy.matrix(np.eye(order,order,1))
    B=cvxpy.matrix(np.vstack((np.zeros((order-1,1)),1)))
    C=(cvxpy.matrix(np.eye(order,order)[:,::-1])*ar).T
    D=cvxpy.matrix([[1]])
    M1=A.T*X
    M2=M1*B
    M=cvxpy.vstack(( \
        cvxpy.hstack((M1*A-X, M2, C.T)), \
        cvxpy.hstack((M2.T, B.T*X*B-H_inf**2, D)), \
        cvxpy.hstack((C, D, cvxpy.matrix([[-1]]))) \
        ))
    constraint1=cvxpy.belongs(-M, cvxpy.semidefinite_cone)
    constraint2=cvxpy.belongs(X,cvxpy.semidefinite_cone)
    p=cvxpy.program(cvxpy.minimize(target),[constraint1, constraint2])
    p.options.update(options)
    quiet = False
    if options.has_key('show_progress'):
        quiet=not options['show_progress']
    p.solve(quiet)
    ntf_ir=np.hstack((1,np.asarray(ar.value.T)[0]))
    return (np.roots(ntf_ir), np.zeros(order), 1.)
