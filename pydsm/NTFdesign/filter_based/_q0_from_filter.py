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

__all__=["q0_from_filter_imp_response"]

def q0_from_filter_imp_response(h_ir, P):
    """
    Computes Q matrix from the output filter impulse response
    
    Parameters
    ----------
    h_ir : array_like
        impulse response of the filter
    P : int
        order of the FIR to be eventually synthesized
        
    Returns
    -------
    q0 : ndarray
        the first row of the matrix Q used in the NTF optimization
        
    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.
    """
    return raw_acorr(h_ir, P)
