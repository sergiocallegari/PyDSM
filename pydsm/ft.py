# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
Fourier transform related routines
==================================

Functions to compute the fft and the dtft.
"""

__all__ = ["fft_centered"]

import numpy as np
import scipy as sp
__import__("scipy.fftpack")
__import__("scipy.integrate")

def fft_centered(x, fs=1):
    """
    Computes a *centered* FFT.
    
    Provides an FFT function where the output vector has the zero frequency 
    at its center. This is more suitable for plotting.
    
    Parameters
    ----------
    x : array_like
        1-D sequence to compute the FFT upon
        
    Returns
    ------
    X : ndarray
        samples of the DFT of the input vector
    ff : ndarray
        vector of frequencies corresponding to the samples in X
        
    Other parameters
    ----------------
    fs : real, optional
        sample frequency for the input vector (defaults to 1)
    """
    N=len(x)
    fs=float(fs)
    if np.mod(N,2)==0:
        k=np.arange(-N/2,N/2) # N even
    else:
        k=np.arange(-(N-1)/2,(N-1)/2+1); # N odd
    ff=k/(N/fs)
    X=sp.fftpack.fft(x)
    X=sp.fftpack.fftshift(X)
    return(ff, X)
