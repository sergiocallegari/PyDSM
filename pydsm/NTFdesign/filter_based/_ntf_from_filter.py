# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
Synthesize a FIR NTF from specs of filter used to remove quantization noise
===========================================================================
"""

import numpy as np
import scipy as sp
__import__("scipy.linalg")
from ._q0_from_filter import q0_from_filter_imp_response
import cvxpy

__all__=["synthesize_ntf_from_filter_ir",
         "synthesize_ntf_from_q0"]

def synthesize_ntf_from_filter_ir(order, h_ir, H_inf=1.5, normalize="auto",
                                 options={}):
    u"""
    Synthesize a FIR NTF based on the ΔΣ modulator output filter.

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
    options : dict, optional
        parameters for the SDP optimizer, see the documentation of `cvxpy`

    Returns
    -------
    ntf : ndarray
        FIR NTF in zpk form
    """
    q0=q0_from_filter_imp_response(h_ir, order)
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
        parameters for the SDP optimizer, see the documentation of `cvxpy`

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
    for opt, val in options:
        p.options[opt]=val
    p.solve()
    ntf_ir=np.hstack((1,np.asarray(ar.value.T)[0]))
    return (np.roots(ntf_ir), np.zeros(order), 1.)
