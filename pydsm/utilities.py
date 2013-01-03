# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
Utility functions for PyDSM
===========================
"""

__all__ = ["is_negligible", "chop", "db", "cplxpair", "mdot"]

import sys
import numpy as np

def is_negligible(x, tol=100*sys.float_info.epsilon):
    """
    Checks if a number is close to zero.
    
    Parameters
    ----------
    x : real
        number to be checked
    tol : real, optional
        absolute tolerance. Defaults to 100 times the system epsilon.
        
    Returns
    -------
    y : bool
        whether the input number is really close to zero or not.
    """
    return abs(np.asarray(x))<tol
    
def chop(x, tol=100*sys.float_info.epsilon):
    """
    Chops to zero the input numbers that are close to zero.

    Parameters
    ----------
    x : real
        number to process
    tol : real, optional
        absolute tolerance. Defaults to 100 times the system epsilon.
        
    Returns
    -------
    y : real
        y is zero if x is close to zero according to the tolerance.
        Alternatively, y is x.
        
    Notes
    -----
    See also `is_negligible`.
    """
    y = np.copy(x)
    y[is_negligible(y, tol)] = 0.
    if np.iscomplexobj(y):
        y = chop(y.real, tol) + 1j*chop(y.imag, tol)
    return y

def db(x, signal_type='voltage', R=1):
    """
    Converts a value to dB a la Matlab
    
    This function tries to replicate the Matlab interface for value to dB
    conversion.    
    
    Parameters
    ----------
    x : real
        value to be converted. Should be positive and non null.
    signal_type : string, optional
        either 'voltage' or 'power'. Defaults to 'voltage'
    R : real, optional
        load resistance value. Defaults to 1.
        
    Returns
    -------
    y : real
        value in dB corresponding to x
        if signal_type is 'power' the result is 10*log10(x). Otherwise (if
        signal_type is 'voltage'), then power is measured over resistor R
        
    Notes
    -----
    The default R value assures that when signal_type is 'voltage' dB defaults
    to the classical 20*log10(x) computation.
    """
    if signal_type == 'power':
        return 10*np.log10(x)
    else:
        return 10*np.log10(np.abs(x)**2/R)
        
def cplxpair(x, tol=100*sys.float_info.epsilon):
    """
    Sorts values in input list by complex pairs.
    
    Parameters
    ----------
    x : array_like of complex
        x is an array of complex values, with the assumption that it contains
        either real values or complex values in conjugate couples.
    tol: real, optional
        absolute tolerance for the recognition of pairs. 
        Defaults to 100 times the system epsilon.

    Returns
    -------
    y : ndarray
        y is an array of complex values, with the same values in x, yet now
        sorted as complex pairs by increasing real part. Real elements in x
        are place after the complex pairs, sorted in increasing order.
        
    Raises
    ------
    ValueError
        'Cannot identify complex pairs.' if there are unpaired complex entries
        in x.
    
    Notes
    -----
    This function is similar to Matlab cplxpair, but not quite.

    """
    x = np.sort_complex(chop(x, tol))
    real_mask = np.isreal(x)
    x_real = x[real_mask]
    x_cplx = x[np.logical_not(real_mask)]
    pos_mask = x_cplx.imag > 0
    x_cplx_pos = x_cplx[pos_mask]
    x_cplx_neg = x_cplx[np.logical_not(pos_mask)]
    x_cplx2 = 1j*np.zeros(2*x_cplx_pos.size)
    for i in range(len(x_cplx_pos)):
        x_cplx2[2*i] = x_cplx_pos[i]
        x_cplx2[2*i+1] = x_cplx_neg[i]
        if abs(x_cplx2[2*i]-np.conj(x_cplx2[2*i+1])) > tol:
            raise ValueError('Cannot identify complex pairs.')
    return np.concatenate((x_cplx2, x_real))

def mdot(*args):
    """
    Dot product taking multiple arguments
    
    Parameters
    ----------
    x1, ..., xn : array_like

    Returns
    -------
    y : ndarray
        y=np.dot(np.dot(np.dot(...np.dot(x1,x2),xn_2)xn_1)xn)    
    """
    return reduce(np.dot, args)
