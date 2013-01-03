# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
Auto and cross correlation
==========================

Functions to compute the auto- and cross- correlation between two vectors.

So far only unnormalized (raw) correlations are supported.
"""

__all__ = ["raw_acorr", "raw_xcorr"]

import numpy as np

def raw_acorr(x, N):
    """
    Computes the raw autocorrelation of a vector up to lag N.
    
    Parameters
    ----------
    x : array_like
        1-D sequence to compute the auto-correlation upon
    N : int
        the maximum (positive) lag of the raw auto-correlation to return.
        
    Returns
    -------
    q : ndarray
        the raw (unnormalized) autocorrelation vector.
        Assuming that m is the length of x
        q(k) = sum_{n=k}^{m-1} x(n) x(n-k) for k = 0 ... N

    Notes
    -----
    The routine does not make any check on the length of x and N. It
    is responsibility of the user to assure that len(x)>=N. In some cases 
    (but only in some cases), zero padding is practiced.
    """
    m=len(x)  
    q=np.asarray([
        np.dot(x[k:m], x[0:m-k])
    for k in xrange(N+1)])    
    return q


def raw_xcorr(x, y, N):
    """
    Computes the raw crosscorrelation between two vectors up to lag N.
    
    Parameters
    ----------
    x : array_like
        first 1-D vector
    y : array_like
        second 1-D vector
    N : int
        the maximum (positive) lag of the raw cross-correlation to return.
        
    Returns
    -------
    q : ndarray
        the raw (unnormalized) crosscorrelation vector.
        Assuming that mx and my are the lengths of x and y
        q(k) = sum_{n=k}^{min(mx-1,my+k-1)} x(n) y(n-k) for k = 0 ... N

    Notes
    -----
    the routine does not make any check on the lengths of x, y and N. It
    is responsibility of the user to assure that N<=len(y). In some cases
    (but only in some cases), zero padding is assumed.
    """
    mx=len(x)
    my=len(y)
    q=np.asarray([
        np.dot(y[k:min(my,mx+k)], x[0:min(my-k,mx)])
    for k in xrange(N+1)])
    return q    
