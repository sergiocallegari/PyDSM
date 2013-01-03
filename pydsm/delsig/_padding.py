# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.
#
# Portions of code ported from the DELSIG toolbox
# Copyright (c) 2009, Richard Schreier

"""
DELSIG style matrix and vector padding routines
===============================================
"""

import numpy as np

__all__=["padl", "padr", "padt", "padb"]

def padl(x, n, val=0):
    """
    Pad a matrix or a vector on the left.
    
    Parameters
    ----------
    x : array_like or matrix_like
        matrix or vector to be padded
    n : int
        (column) length to pad x to
    val : number, optional
        fill value for padding. Defaults to 0.
        
    Returns
    -------
    y : ndarray
        matrix or vector x, padded to the left to have n entries (if vector)
        or n columns (if matrix). Padding value is val.
        
    Notes
    -----
    The empty matrix is assumed to have 1 empty column.
    """
    if x.ndim == 1:
        return np.hstack((val*np.ones(n-x.shape[0]), x))
    else:
        return np.hstack((val*np.ones((x.shape[0], n-x.shape[1])), x))

def padr(x, n, val=0):
    """Pad a matrix or a vector x on the right.
      
    Parameters
    ----------
    x : array_like or matrix_like
        matrix or vector to be padded
    n : int
        (column) length to pad x to
    val : number, optional
        fill value for padding. Defaults to 0.
        
    Returns
    -------
    y : ndarray
        matrix or vector x, padded to the right to have n entries (if vector)
        or n columns (if matrix). Padding value is val.
        
    Notes
    -----
    The empty matrix is assumed to have 1 empty column.
    """
    if x.ndim == 1:
        return np.hstack((x, (val*np.ones(n-x.shape[0]))))
    else:
        return np.hstack((x, (val*np.ones((x.shape[0], n-x.shape[1])))))
        
def padt(x, n, val=0):
    """Pad a matrix x on the top.
    
    Parameters
    ----------
    x : matrix_like
        matrix to be padded
    n : int
        row height to pad x to
    val : number, optional
        fill value for padding. Defaults to 0.
        
    Returns
    -------
    y : ndarray
        matrix x, padded to the top to have n rows. Padding value is val.
        
    Notes
    -----
    The empty matrix is assumed to have 1 empty column.
    """
    return np.vstack((val*np.ones((n-x.shape[0], x.shape[1])), x))

def padb(x, n, val=0):
    """Pad a matrix x on the bottom.
    
    Parameters
    ----------
    x : matrix_like
        matrix to be padded
    n : int
        row height to pad x to
    val : number, optional
        fill value for padding. Defaults to 0.
        
    Returns
    -------
    y : ndarray
        matrix x, padded to the bottom to have n rows. Padding value is val.
        
    Notes
    -----
    The empty matrix is assumed to have 1 empty column.
    """
    return np.vstack((x, val*np.ones((n-x.shape[0], x.shape[1]))))
