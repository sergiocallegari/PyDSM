# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

# This file is part of PyDSM.

# PyDSM is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# PyDSM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with PyDSM.  If not, see <https://www.gnu.org/licenses/>.

# This file includes code ported from the DELSIG Matlab toolbox
# (see https://www.mathworks.com/matlabcentral/fileexchange/19)
# covered by the following copyright and permission notice
#
# Copyright (c) 2009 Richard Schreier
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the distribution
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""
Collection of DELSIG style plotting routines
============================================
"""

import numpy as np
import matplotlib.pyplot as plt
from ..relab import cplxpair
from scipy.signal import tf2zpk

__all__ = ["plotPZ"]


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

    Other Parameters
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
        H = tf2zpk(*H)
    z = cplxpair(H[0])
    p = cplxpair(H[1])

    plt.grid(True)

    # Plot x and o for poles and zeros, respectively
    plt.plot(p.real, p.imag, linestyle='None', **pole_fmt)
    if len(z) > 0:
        plt.plot(z.real, z.imag, linestyle='None', **zero_fmt)

    # Draw unit circle, real axis and imag axis
    circle = np.exp(2j*np.pi*np.linspace(0, 1, 100))
    plt.plot(circle.real, circle.imag, 'k')
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
                 horizontalalignment='left',
                 verticalalignment='center')
        y = y - 0.1
        for i in range(len(pp)):
            if pp[i].imag == 0:
                str_p = '$%+.4f$' % pp[i].real
            else:
                str_p = r'$%+.4f\pm \mathrm{j}%.4f$' % (pp[i].real, pp[i].imag)
            plt.text(-0.9, y, str_p,
                     horizontalalignment='left',
                     verticalalignment='center')
            y = y - 0.1
        if len(z) > 0:
            zz = z[z.imag >= 0]
            y = 0.05*(len(zz)+1)
            str_z = 'Zeros:'
            plt.text(0, y, str_z,
                     horizontalalignment='left',
                     verticalalignment='center')
            y = y - 0.1
            for i in range(len(zz)):
                if zz[i].imag == 0:
                    str_z = '$%+.4f$' % zz[i].real
                else:
                    str_z = (r'$%+.4f\pm \mathrm{j}%.4f$' %
                             (zz[i].real, zz[i].imag))
                plt.text(0, y, str_z,
                         horizontalalignment='left',
                         verticalalignment='center')
                y = y - 0.1

    plt.ylabel('Imag')
    plt.xlabel('Real')
