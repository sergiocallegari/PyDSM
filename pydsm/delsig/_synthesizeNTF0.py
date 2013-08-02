# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.
#
# Portions of code ported from the DELSIG toolbox
# Copyright (c) 2009, Richard Schreier

import numpy as np
from warnings import warn
from ..errors import PyDsmWarning, PyDsmError
from ._tf import evalTF
from ..utilities import cplxpair
from ._ds import ds_optzeros

def synthesizeNTF0(order, osr, opt, H_inf, f0):
    # Determine the zeros.
    if f0 != 0:
        # Bandpass design-- halve the order temporarily.
        order = order/2
        dw = np.pi/(2*osr)
    else:
        dw = np.pi/osr

    if opt.ndim == 0:
        # opt is a number
        if opt == 0:
            z = np.zeros(order)
        else:
            z = dw*ds_optzeros(order, opt)
        if z.size == 0:
            raise PyDsmError('Cannot synthesize NTF zeros')
        if f0 != 0:
            # Bandpass design-- shift and replicate the zeros.
            order = order*2
            z = z + 2*np.pi*f0
            z = np.vstack((z,-z)).transpose().flatten()
        z = np.exp(1j*z)
    else:
        z = opt

    p = np.zeros(order)
    k = 1
    itn_limit = 100
    fprev = 0

    if f0 == 0:
        # Lowpass design
        HinfLimit = 2**order
        # !!! The limit is actually lower for opt=1 and low OSR
        if H_inf >= HinfLimit:
            warn('Unable to achieve specified H_inf.\n'
                'Setting all NTF poles to zero.',
                PyDsmWarning)
            p = np.zeros(order)
        else:
            x = 0.3**(order-1)   # starting guess
            for itn in xrange(1, itn_limit+1):
                me2 = -0.5*(x**(2./order))
                w = (2*np.arange(1,order+1)+1)*np.pi/order
                mb2 = 1+me2*np.exp(1j*w)
                p = mb2 - np.sqrt(mb2**2-1)
                # Reflect poles to be inside the unit circle
                out = abs(p)>1
                p[out] = 1/p[out]
                # The following is not exactly what delsig does.
                # We do not have an identical cplxpair
                p = cplxpair(p)
                f = np.real(evalTF((z, p, k), -1))-H_inf
                if itn == 1:
                    delta_x = -f/100
                else:
                    delta_x = -f*delta_x/(f-fprev)
                xplus = x+delta_x
                if xplus > 0:
                    x = xplus
                else:
                    x = x*0.1
                fprev = f
                if abs(f)<1e-10 or abs(delta_x)<1e-10:
                    break
                if x > 1e6:
                    warn('Unable to achieve specified Hinf.\n'
                         'Setting all NTF poles to zero.', PyDsmWarning)
                    p = np.zeros(order)
                    break
                if itn == itn_limit:
                    warn('Danger! Iteration limit exceeded.', PyDsmWarning)
    else:
        # Bandpass design
        x = 0.3**(order/2-1)   # starting guess (not very good for f0~0)
        if f0 > 0.25:
            z_inf = 1.
        else:
            z_inf = -1.
        c2pif0 = np.cos(2*np.pi*f0)
        for itn in xrange(1, itn_limit+1):
            e2 = 0.5*x**(2./order)
            w = (2*np.arange(order)+1)*np.pi/order
            mb2 = c2pif0 + e2*np.exp(1j*w)
            p = mb2 - np.sqrt(mb2**2-1)
            # Reflect poles to be inside the unit circle
            out = abs(p)>1
            p[out] = 1/p[out]
            # The following is not exactly what delsig does.
            p = cplxpair(p)
            f = np.real(evalTF((z, p, k), z_inf))-H_inf
            if itn == 1:
                delta_x = -f/100
            else:
                delta_x = -f*delta_x/(f-fprev)
            xplus = x+delta_x
            if xplus > 0:
                x = xplus
            else:
                x = x*0.1
            fprev = f
            if abs(f)<1e-10 or abs(delta_x)<1e-10:
                break
            if x > 1e6:
                warn('Unable to achieve specified Hinf.\n'
                    'Setting all NTF poles to zero.', PyDsmWarning)
                p = np.zeros(order)
                break
            if itn == itn_limit:
                warn('Danger! Iteration limit exceeded.', PyDsmWarning)

    z = cplxpair(z)
    return z, p, k
