# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

# Portions of code ported from the DELSIG toolbox
# Copyright (c) 2009, Richard Schreier

"""
Transfer function management a la DELSIG.
=========================================
"""

import numpy as np

__all__=["evalTF", "evalRPoly", "rmsGain"]

def evalTF(tf, x):
    """
    Evaluates a transfer function.
    
    Parameters
    ----------
    tf : tuple
        transfer function in zpk or nd form
    x : complex or array_like of complex
        value or vector of values where the tf is to be evaluated
    
    Returns
    -------
    y : ndarray
        value(s) of transfer function at the given complex values.
        
    Note
    ----
    x corresponds to 's' or 'z' in CT or DT transfer functions respectively.
    Thus it should be 1j*omega or exp(1j*omega*T).
    """
    if len(tf) == 3:
        return evalRPoly(tf[0], x, tf[2])/evalRPoly(tf[1], x, 1)
    elif len(tf) == 2:
        return np.polyval(tf[0], x)/np.polyval(tf[1],x)
        
def evalRPoly(roots, x, k=1):
    """
    Compute the value of a polynomial that is given in terms of its roots.
    
    Parameters
    ----------
    roots : array_like
        roots of polynomial
    x : complex or array_like of complex
        complex value or vector of complex values where tf is to be
        evaluated
    k -> real, optional
        gain, defaults to 1.
        
    Returns
    -------
    y : ndarray
        value of polynomial at the given complex values.        
    """
    y = k*np.ones_like(x)
    roots = np.asarray(roots)
    # Remove roots at infinity
    roots = roots[np.logical_not(np.isinf(roots))]
    for i in xrange(roots.size):
        y = y*(x-roots[i])
    return y
    
def rmsGain(H, f1, f2, N=100):
    """
    Compute the root mean-square gain of a DT transfer function.
    
    Parameters
    ----------
    H : tuple
        transfer function either in (z,p,k) or (n,d) form
    f1 : real
        lower bound of frequency band on which the transfer function
        is evaluated
    f2 : real 
        upper bound of frequency band on which the transfer function
        is evaluated
    N : int, optional
        number of points where the transfer function is evaluated in the
        interval. Defaults to 100.
            
    Returns
    -------
    rms : real
        rms value of the discrete time transfer function.
    
    Notes
    -----
    The discrete-time transer function H is evaluated in the frequency band
    (f1,f2).  Spanning of the bandwidth is linear. Frequencies are normalized
    in the [0,0.5] interval.
    """
    w = np.linspace(2*np.pi*f1, 2*np.pi*f2, N)
    return np.linalg.norm( evalTF(H, np.exp(1j*w))) / np.sqrt(N)
