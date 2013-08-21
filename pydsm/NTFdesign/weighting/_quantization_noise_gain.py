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

__all__ = ["quantization_weighted_noise_gain"]

def quantization_weighted_noise_gain(NTF, w):
    r"""
    Computes the quantization noise power gain

    Parameters
    ----------
    NTF : tuple
        NTF definition in zpk or nd form
    w : callable with argument f in [0,1/2]
        noise weighting function

    Returns
    -------
    a : real
        noise power gain

    Notes
    -----
    The computation is practiced as

    .. math::
        \int_{f=0}^{\frac{1}{2}}
        \left|\mathit{NTF}
        \left(\mathrm{e}^{\mathrm{i} 2\pi f}\right)\right|^2 w(f) df
    """
    return sp.integrate.quad(lambda f: \
        np.abs(evalTF(NTF,np.exp(2j*np.pi*f)))**2 * w(f), 0, 0.5)[0]