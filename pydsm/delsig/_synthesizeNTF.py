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
# along with PyDSM.  If not, see <http://www.gnu.org/licenses/>.

# This file includes code ported from the DELSIG Matlab toolbox
# (see http://www.mathworks.com/matlabcentral/fileexchange/19)
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
Entry point for DELSIG-like Delta-Sigma NTF synthesis function
==============================================================
"""

import numpy as np
from warnings import warn
from ..exceptions import PyDsmApproximationWarning
from ._synthesizeNTF0 import synthesizeNTF0
from ._synthesizeNTF1 import synthesizeNTF1

__all__ = ["optimize_NTF", "synthesizeNTF"]

optimize_NTF = True

# filterwarnings("always", ".*", Warning,
#                "pydsm.synthesis._synthesizeNTF[01]?")


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
        * 3 -> optimized zeros (with optimizer)
        * 4 -> same as 3, but with at least one zero at band-center
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
        'Frequency f0 must be less than 0.5' if f0 is out of range

        'Order must be even for a bandpass modulator' if the order is
        incompatible with the modulator type.

        'The opt vector must be of length xxx' if opt is used to explicitly
        pass the NTF zeros and these are in the wrong number.

    RuntimeError
        'Cannot synthesize NTF zeros' if the synthesis fails for some
        unspecified reason.

    Warns
    -----
    PyDsmApproximationWarning
        'Creating a lowpass ntf.' if the center frequency is different
        from zero, but so low that a low pass modulator must be designed.

        'Unable to achieve specified H_inf ...' if the desired H_inf
        cannot be achieved.

        'Iteration limit exceeded' if the routine converges too
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
        raise ValueError('Frequency f0 must be less than 0.5')
    if f0 != 0 and f0 < 0.25/osr:
        warn('Creating a lowpass ntf.', PyDsmApproximationWarning)
        f0 = 0
    if f0 != 0 and order % 2 != 0:
        raise ValueError('Order must be even for a bandpass modulator')
    opt = np.asarray(opt)
    if opt.ndim > 1 or (opt.ndim == 1 and opt.size != order):
        raise ValueError('The opt vector must be of length %d' % order)

    if not optimize_NTF:
        ntf = synthesizeNTF0(order, osr, opt, H_inf, f0)
    else:
        ntf = synthesizeNTF1(order, osr, opt, H_inf, f0)
    return ntf
