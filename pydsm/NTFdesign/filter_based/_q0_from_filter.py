# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
Functions to obtain the matrix used in the NTF optimization
===========================================================

The matrix is symmetric Toeplitz, thus actually described by its first
row only, which is what the routines here return.
"""

from ...correlations import raw_acorr
from ..weighting import q0_from_noise_weighting
from ...delsig import evalTF
import numpy as np
from warnings import warn
from ...errors import PyDsmError, PyDsmWarning

__all__=["q0_from_filter",
         "q0_from_filter_imp_response",
         "q0_from_filter_mag_response"]


def q0_from_filter_imp_response(P, h_ir):
    """
    Computes Q matrix from the output filter impulse response

    Parameters
    ----------
    P : int
        order of the FIR to be eventually synthesized
    h_ir : array_like
        impulse response of the filter

    Returns
    -------
    q0 : ndarray
        the first row of the matrix Q used in the NTF optimization

    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.
    """
    warn('Deprecated function q0_from_filter_imp_response.\n'
        'Will be removed shortly.\n'
        'Use q0_from_filter instead.', PyDsmWarning)
    return q0_from_filter(P, h_ir, F_type='imp')


def q0_from_filter_mag_response(P, h_mag,\
    integrator_params={'epsabs':1E-14, 'epsrel':1E-9}):
    """
    Computes Q matrix from the output filter magnitude response

    Parameters
    ----------
    P : int
        order of the FIR to be eventually synthesized
    h_mag : callable
        function of f representing the filter magnitude response
        f is normalized between 0 and 0.5

    Returns
    -------
    q0 : ndarray
        the first row of the matrix Q used in the NTF optimization

    Other parameters
    ----------------
    integrator_params : dict, optional
        the controlling parameters for the numerical integrator
        (see `scipy.integrate.quad`)

    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.
    """
    warn('Deprecated function q0_from_filter_mag_response.\n'
        'Will be removed shortly.\n'
        'Use q0_from_filter instead.', PyDsmWarning)
    return q0_from_filter(P, h_mag, F_type='mag',
                          integrator_params=integrator_params)


def q0_from_filter(P, F, F_type,
                   integrator_params={'epsabs':1E-14, 'epsrel':1E-9}):
    """
    Computes Q matrix from the modulator output filter

    Parameters
    ----------
    P : int
        order of the FIR to be eventually synthesized
    F : tuple or array_like or callable
        output filter description.
        This is given by a zpk or nd form if F_type is 'zpk' or 'nd'.
        It is a magnitude response (function of f, with f in [0,1/2]) if
        F_type is 'mag'.
        It is an impulse response if F_type is 'imp'.
    F_type : str
        string indicating how F is expressed. Can be 'zpk', 'nd', 'mag' or
        'imp'

    Returns
    -------
    q0 : ndarray
        the first row of the matrix Q used in the NTF optimization

    Other parameters
    ----------------
    integrator_params : dict, optional
        the controlling parameters for the numerical integrator
        (see `scipy.integrate.quad`). Unused if F_type is 'imp'.

    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.
    """
    if F_type=='zpk' or F_type=='nd':
        w = lambda f: np.abs(evalTF(F,np.exp(2j*np.pi*f)))**2
        q0 = q0_from_noise_weighting(P, w, integrator_params)
    elif F_type=='imp':
        q0 = raw_acorr(F, P)
    elif F_type=='mag':
        w = lambda f: F(f)**2
        q0 = q0_from_noise_weighting(P, w, integrator_params)
    else:
        raise PyDsmError("Incorrect filter type specification")
    return q0
