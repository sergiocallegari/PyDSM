# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
Synthesize a FIR NTF from specs of filter used to remove quantization noise
===========================================================================
"""

from ._q0_from_filter import (q0_from_filter_imp_response,
                              q0_from_filter_mag_response)
from ..weighting import synthesize_ntf_from_q0
from warnings import warn
from ...errors import PyDsmWarning

__all__=["synthesize_ntf_from_filter_imp",
         "synthesize_ntf_from_filter_mag",
         "synthesize_ntf_from_filter_ir"]

def synthesize_ntf_from_filter_ir(order, h_ir, H_inf=1.5, normalize="auto",
                                  options={}):
    warn('Deprecated function synthesize_ntf_from_filter_ir.\n'
        'Will be removed shortly.\n'
        'Use synthesize_ntf_from_filter_imp instead.', PyDsmWarning)
    return synthesize_ntf_from_filter_imp(order, h_ir, H_inf, normalize,
                                          options)

def synthesize_ntf_from_filter_imp(order, h_ir, H_inf=1.5, normalize="auto",
                                   options={}):
    u"""
    Synthesize a FIR NTF based on the ΔΣM output filter impulse response.

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
        parameters for the SDP optimizer, see the documentation of `cvxpy`.
        This includes 'show_progress' (default True).

    Returns
    -------
    ntf : ndarray
        FIR NTF in zpk form
    """
    q0=q0_from_filter_imp_response(order, h_ir)
    return synthesize_ntf_from_q0(q0, H_inf, normalize, options)

def synthesize_ntf_from_filter_mag(order, h_mag, H_inf=1.5, normalize="auto",
                                  options={}):
    u"""
    Synthesize a FIR NTF based on the ΔΣM output filter magnitude response.

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
    options : dict, optional
        parameters for the SDP optimizer, see the documentation of `cvxpy`.
        This includes 'show_progress' (default True).

    Returns
    -------
    ntf : ndarray
        FIR NTF in zpk form

    Notes
    -----
    The computation of the NTF from the output filter magnitude response
    involves computing an integral on the magnitude response. To control the
    integration parameters, do not use this function. Rather, first compute a
    vector q0 with `q0_from_filter_mag_response` (which lets the integrator
    params be specified), then use `synthesize_ntf_from_q0`.
    """
    q0=q0_from_filter_mag_response(order, h_mag)
    return synthesize_ntf_from_q0(q0, H_inf, normalize, options)
