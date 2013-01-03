# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
Computation of the noise power gain through the NTF and the output filter
=========================================================================
"""

import numpy as np
import scipy as sp
__import__("scipy.integrate")
from ...ir import impulse_response
from ...delsig import evalTF

__all__ = ["quantization_noise_gain", "quantization_noise_gain_ir"]

def quantization_noise_gain(NTF, H):
    """
    Computes the quantization noise power gain
    
    Parameters
    ----------
        NTF : tuple
            NTF definition in zpk or nd form
        H : tuple
            output filter definition in zpk or nd form
            
    Returns
    -------
        a : real
            noise power gain
    
    Notes
    -----
    The computation is practiced as
        int_{w=0}^{inf} |H(jw)|^2 |NTF(jw)|^2 dw
    """
    def fprod(h1, h2, f):
        return np.abs(evalTF(h1,np.exp(1j*2*np.pi*f)))**2* \
            np.abs(evalTF(h2,np.exp(1j*2*np.pi*f)))**2
    return sp.integrate.quad(lambda f: fprod(NTF, H, f), 0, 0.5)[0]
    
def quantization_noise_gain_ir(NTF, H):
    """
    Computes the quantization noise power gain
    
    Parameters
    ----------
        NTF : tuple
            NTF definition in zpk or nd form
        H : tuple
            output filter definition in zpk or nd form
            
    Returns
    -------
        a : real
            noise power gain
    
    Notes
    -----
    The computation is practiced as the sum of the squared entries
    in the impulse response of the cascaded filter NTF*H
    """    
    h1_ir=impulse_response(NTF)
    h2_ir=impulse_response(H)
    conv = sp.signal.convolve(h1_ir, h2_ir)
    return np.sum(conv**2)/2
