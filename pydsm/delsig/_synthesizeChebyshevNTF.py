# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.
#
# Portions of code ported from the DELSIG toolbox
# Copyright (c) 2009, Richard Schreier

"""
DELSIG-like synthesizeChebyshevNTF
==================================
"""

import numpy as np
from scipy.signal import cheby2
from ._ds import ds_f1f2
from ..utilities import cplxpair

__all__=["synthesizeChebyshevNTF"]


def synthesizeChebyshevNTF(order=3, osr=64, opt=0, H_inf=1.5, f0=0.0):
    """
    Synthesize the NTF for a DS modulator by a Chebyshev function.

    This function is an alternative to synthesizeNTF which assumes
    that the magnitude of the denominator of the NTF is approximately constant
    in the passband. When the OSR or H_inf are low, this assumption breaks
    down and synthesizeNTF yields a non-optimal NTF. The synthesizeChebyshevNTF
    creates non-optimal NTFs, but fares better than synthesizeNTF in the
    aforementioned circumstances. The NTF is obtained from a type-2 highpass
    Chebyshev function.

    Parameters
    ----------
    order : int, optional
        the order of the modulator, defaults to 3
    osr : float, optional
        the oversamping ratio (based on the actual signal bandwidth)
    opt : int, optional
        ignored, present for signature compatibility with synthesizeNTF
    H_inf : real, optional
        max allowed peak value of the NTF. Defaults to 1.5
    f0 : real, optional
        center frequency, defaults to 0.
        1 corresponds to the sampling frequency, so that 0.5 is the
        maximum value. Value 0 specifies an LP modulator.

    Returns
    -------
    ntf : tuple
        noise transfer function in zpk form.

    Raises
    ------
    ValueError
        'Order must be even for a bandpass modulator.' if the order is
        incompatible with the modulator type.

    Notes
    -----
    Parameter opt is ignored.
    """
    if f0 != 0 and order % 2 != 0:
        raise ValueError('Order must be even for a bandpass modulator.')
    else:
        f1, f2 = ds_f1f2(osr,f0)
        f1f2 = (f1, f2)

    # Iteratively solve for the attenuation spec (x) which yields the
    # desired H_inf
    x_min = 0.
    x_max = 300.
    dx_max = 10.
    ftol = 1e-6
    xtol = 1e-6
    # Initial guess for x
    x = 60.
    itn_limit=10
    converged = False
    f_p = 0.

    # Design is based on a Chebyshev type II filter, designed with the cheby2
    # function. Recall that a Chebyshev type II filter is monothonic in the
    # pass-band and equiripple in the stop-band.
    # In cheby2 parameters are:
    #   N : order
    #   Rs : stopband attenuation
    #        this is measured from the peak value in the pass-band
    #        to the peak(s) of the stop-band, in dB
    #   Wn : stopband edge frequency (normalized so that 1 is fs/2 for DT
    #        filters). Must be a couple of frequencies for stop band filters
    #   btype : filter type (e.g. 'low', 'high', 'stop')
    #   analog : non zero for a CT filter (defaults to DT)
    #   output : 'ba' or 'zpk'
    # The filter that is returned has peak gain at 1 in the pass-band

    for itn in xrange(itn_limit):
        if f0 == 0:
            z, p, k = cheby2(order, x, 1./osr, 'high', output='zpk')
        else:
            z, p, k = cheby2(order/2, x, 2*f1f2, 'stop',
                                       output='zpk')
        f = 1./k - H_inf
        # print (x, f)
        if f>0:
            # x is too big
            x_max = x
        else:
            # x is too small
            x_min = x
        if itn==0:
            # First iteration
            dx = -dx_max * np.sign(f)
        else:
            # Normal iteration
            df = f - f_p
            if np.abs(df) < ftol:
                converged = True
                break
            dx = -f * dx/df;
        if converged:
            break
        x_p = x
        f_p = f
        x = np.max((x_min, np.min((x+dx, x_max))))
        dx = x - x_p
        if np.abs(dx) < xtol:
            break

    z=cplxpair(z)
    p=cplxpair(p)
    return (z, p, 1)
