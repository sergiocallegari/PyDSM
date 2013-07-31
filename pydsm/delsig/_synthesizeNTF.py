# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.
#
# Portions of code ported from the DELSIG toolbox
# Copyright (c) 2009, Richard Schreier

"""
Entry point for DELSIG-like Delta-Sigma NTF synthesis function
==============================================================
"""

import numpy as np
from warnings import warn, filterwarnings
from ..errors import PyDsmWarning
from ._synthesizeNTF0 import synthesizeNTF0
from ._synthesizeNTF1 import synthesizeNTF1

__all__=["optimize_NTF", "synthesizeNTF"]

optimize_NTF = True

#filterwarnings("always", ".*", Warning,
#               "pydsm.synthesis._synthesizeNTF[01]?")

def synthesizeNTF(order=3, osr=64, opt=0, H_inf=1.5, f0=0.0):
    """
    Synthesizes an NTF for a DS modulator by Schreier's approach.

    Parameters
    ----------
    order : int, optional
        the order of the modulator, defaults to 3
    osr : float, optional
        the oversamping ratio (based on the actual signal bandwidth)
    opt : int or list of floats, optional
        flag for optimized zeros, defaults to 0

        * 0 -> not optimized,
        * 1 -> optimized,
        * 2 -> optimized with at least one zero at band-center,
        * 3 -> optimized zeros (Requires MATLAB6 and Optimization Toolbox)
        * [z] -> zero locations in complex form

    H_inf : real, optional
        max allowed peak value of the NTF. Defaults to 1.5
    f0 : real, optional
        center frequency for BP modulators, or 0 for LP modulators.
        Defaults to 0.
        1 corresponds to the sampling frequency, so that 0.5 is the
        maximum value. Value 0 specifies an LP modulator.

    Returns
    -------
    ntf : tuple
        noise transfer function in zpk form.

    Raises
    ------
    ValueError
        'Error. f0 must be less than 0.5' if f0 is out of range

        'Order must be even for a bandpass modulator.' if the order is
        incompatible with the modulator type.

        'The opt vector must be of length xxx' if opt is used to explicitly
        pass the NTF zeros and these are in the wrong number.

    Warns
    -----
    PyDsmWarning
        'Creating a lowpass ntf.' if the center frequency is different
        from zero, but so low that a low pass modulator must be designed.

        'Unable to achieve specified H_inf ...' if the desired H_inf
        cannot be achieved.

        'Danger! Iteration limit exceeded' if the routine converges too
        slowly.

    Notes
    -----
    This is actually a wrapper function which calls the appropriate version
    of synthesizeNTF, based on the control flag `optimize_NTF` which
    determines whether to use optimization tools.

    Parameter H_inf is used to enforce the Lee stability criterion.

    If osr or H_inf are low, then the NTF is non optimal. Use
    synthesizeChebyshevNTF instead.
    """
    if f0 > 0.5:
        raise ValueError('Error. f0 must be less than 0.5.')
    if f0 != 0 and f0 < 0.25/osr:
        warn('Creating a lowpass ntf.', PyDsmWarning)
        f0 = 0
    if f0 != 0 and order % 2 != 0:
        raise ValueError('Order must be even for a bandpass modulator.')
    opt = np.asarray(opt)
    if opt.ndim > 1 or (opt.ndim == 1 and opt.size != order):
        raise ValueError('The opt vector must be of length %d.' % order)

    if optimize_NTF == False:
        ntf = synthesizeNTF0(order, osr, opt, H_inf, f0)
    else:
        ntf = synthesizeNTF1(order, osr, opt, H_inf, f0)
    return ntf
