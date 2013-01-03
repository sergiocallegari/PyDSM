# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.
#
# Portions of code ported from the DELSIG toolbox
# Copyright (c) 2009, Richard Schreier

"""
Collection of DELSIG style plotting routines
============================================
"""

import numpy as np
import matplotlib.pyplot as plt
from ..utilities import cplxpair
from scipy.signal import tf2zpk

__all__=["plotPZ"]

def plotPZ(H, color='b', markersize=5, showlist=False):
    """
    Plots the poles and zeros of a transfer function.
    
    Parameters
    ----------
    H : tuple
        transfer function in pzk or nd form
    showlist : bool
        if showlist is true, a list of the poles and zeros is superimposed
        onto the plot.
        
    Other parameters
    ----------------
    color : string or list of strings, optional
        color or colors to plot the poles and the zeros (defaults to 'b'
        meaning black)
    markersize : real, optional
        size of the markers used to represent the poles and the zeros
        (defaults to 5)
        
    Notes
    -----
    See `matplotlib` for info about color codes.        
    
    """

    pole_fmt = {'marker': 'x', 'markersize': markersize}
    zero_fmt = {'marker': 'o', 'markersize': markersize}

    if isinstance(color, list):
        pole_fmt['color'] = color[0]
        zero_fmt['color'] = color[1]
    else:
        pole_fmt['color'] = color
        zero_fmt['color'] = color

    if len(H) == 2:
        H = tf2zpk(**H)
    z = cplxpair(H[0])
    p = cplxpair(H[1])
 
    hold_status = plt.ishold()
    plt.grid(True)

    # Plot x and o for poles and zeros, respectively    
    plt.plot(p.real, p.imag, linestyle='None', **pole_fmt)
    plt.hold(True)
    if len(z) > 0:
        plt.plot(z.real, z.imag, linestyle='None', **zero_fmt)
        
    # Draw unit circle, real axis and imag axis
    circle = np.exp(2j*np.pi*np.linspace(0, 1, 100))
    plt.plot(circle.real, circle.imag)
    plt.axis('equal')
    limits = plt.axis()
    plt.plot([0, 0], limits[1:3], 'k:')
    plt.plot(limits[0:2], [0, 0], 'k:')

    if showlist:
        # List the poles and zeros
        pp = p[p.imag >= 0]
        y = 0.05*(len(pp)+1)
        str_p = 'Poles:'
        plt.text(-0.9, y, str_p,
                 horizontalalignment = 'left',
                 verticalalignment = 'center')
        y = y - 0.1
        for i in xrange(len(pp)):
            if pp[i].imag == 0:
                str_p = '%+.4f' % pp[i].real
            else:
                str_p = '%+.4f+/-j%.4f' %  (pp[i].real, pp[i].imag)
            plt.text(-0.9, y, str_p,
                     horizontalalignment = 'left',
                     verticalalignment = 'center')
            y = y - 0.1
        if len(z) > 0:
            zz = z[z.imag >= 0]
            y = 0.05*(len(zz)+1)
            str_z = 'Zeros:'
            plt.text(0, y, str_z,
                     horizontalalignment = 'left',
                     verticalalignment = 'center')
            y = y - 0.1
            for i in xrange(len(zz)):
                if zz[i].imag == 0:
                    str_z = '%+.4f' % zz[i].real
                else:
                    str_z = '%+.4f+/-j%.4f' % (zz[i].real, zz[i].imag)
                plt.text(0, y, str_z, 
                         horizontalalignment = 'left',
                         verticalalignment = 'center')
                y = y - 0.1
     
    if not hold_status:
        plt.hold(False)
