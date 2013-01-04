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
from _q0_from_filter import q0_from_filter_imp_response
import cvxpy

__all__=["synthesize_ntf_from_filter_ir",
         "synthesize_ntf_from_q0"]

def synthesize_ntf_from_filter_ir(h_ir, P, gamma=1.5, mult="auto",
                                 options={}):
    """
    Synthesize a FIR NTF based on the impulse response of the filter
    in charge of removing the quantization noise

    Parameters
    ----------
    h_ir : array_like of reals
        filter impulse response
    P : int
        Delta sigma modulator order
    gamma : real, optional
        Max peak NTF gain, defaults to 1.5, used to enforce the Lee criterion
    mult : string or real, optional
        Normalization to apply to the quadratic form used in the NTF
        selection. Defaults to 'auto' which means setting the top left entry
        in the matrix Q defining the quadratic form to 1.
    options : dict, optional
        parameters for the SDP optimizer, see the documentation of `cvxpy`

    Returns
    -------
    ntf : ndarray
        coefficients of the FIR NTF
    """
    q0=q0_from_filter_imp_response(h_ir, P)
    return synthesize_ntf_from_q0(q0, gamma, mult, options)

def synthesize_ntf_from_q0(q0, gamma=1.5, mult="auto",
                          options={}):
    """
    Synthesize a FIR NTF based a quadratic form defining the quantization
    noise gain

    Parameters
    ----------
    q0 : array_like
        first row of the Toeplitz symmetric matrix defining the quadratic form
    gamma : real, optional
        Max peak NTF gain, defaults to 1.5, used to enforce the Lee criterion
    mult : string or real, optional
        Normalization to apply to the quadratic form used in the NTF
        selection. Defaults to 'auto' which means setting the top left entry
        in the matrix Q defining the quadratic form to 1.
    options : dict, optional
        parameters for the SDP optimizer, see the documentation of `cvxpy`

    Returns
    -------
    ntf : ndarray
        coefficients of the FIR NTF
    """
    P=q0.shape[0]-1
    if mult=='auto':
        q0=q0/q0[0]
    elif mult!=None:
        q0=q0*mult
    Q=cvxpy.matrix(sp.linalg.toeplitz(q0[0:-1]))
    L=cvxpy.matrix(2*q0[1:])
    ar=cvxpy.variable(P,1, name='ar')
    target=cvxpy.quad_form(ar,Q) + L*ar
    X=cvxpy.variable(P, P, structure='symmetric', name='X')
    A=cvxpy.matrix(np.eye(P,P,1))
    B=cvxpy.matrix(np.vstack((np.zeros((P-1,1)),1)))
    C=(cvxpy.matrix(np.eye(P,P)[:,::-1])*ar).T
    D=cvxpy.matrix([[1]])
    # If X is symmetric, then M is symmetric ?
    M1=A.T*X
    M2=M1*B
    M=cvxpy.vstack(( \
        cvxpy.hstack((M1*A-X, M2, C.T)), \
        cvxpy.hstack((M2.T, B.T*X*B-gamma**2, D)), \
        cvxpy.hstack((C, D, cvxpy.matrix([[-1]]))) \
        ))
    constraint1=cvxpy.belongs(-M, cvxpy.semidefinite_cone)
    constraint2=cvxpy.belongs(X,cvxpy.semidefinite_cone)
    p=cvxpy.program(cvxpy.minimize(target),[constraint1, constraint2])
    for opt, val in options:
        p.options[opt]=val
    p.solve()
    ntf_ir=np.hstack((1,np.asarray(ar.value.T)[0]))
    return ntf_ir
