# -*- coding: utf-8 -*-

# Copyright (c) 2013, Sergio Callegari
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
CLANS NTF design method
=======================

CLANS is the Closed-Loop Analysis of Noise Shaping Coders method
proposed by J. G. Kenney and L. R. Carley in
'Design of Multibit Noise Shaping Data Converters'
Analog Integrated Circuits and Signal Processing Vol. 3 pp 259-272 (1993)
"""

import numpy as np
from scipy.optimize import minimize
from ._synthesizeNTF import synthesizeNTF
from ..relab import cplxpair
from ._tf import evalTF
from ..ir import impulse_response
from ._dsclansNTF import dsclansNTF

__all__ = ["clans"]


def clans(order=4, osr=64, nq=5, rmax=0.95, opt=0, **options):
    u"""Synthesize the NTF for a ΔΣM w/ multibit quantizer by the CLANS method.

    This function is based on the CLANS method (Closed-Loop Analysis of Noise
    Shaping Coders) by J. G. Kenney and L. R. Carley, described in
    'Design of Multibit Noise Shaping Data Converters'
    Analog Integrated Circuits and Signal Processing Vol. 3 pp 259-272 (1993).

    The method is only suitable for low pass modulators.

    Parameters
    ----------
    order : int, optional
        the order of the modulator, defaults to 4
    osr : float, optional
        the oversamping ratio (based on the actual signal bandwidth).
        Defaults to 64
    nq : int, optional
        the number of quantization levels, defaults to 5
    rmax : real, optional
        the maximum radius of the NTF poles
    opt : int, optional
        opt parameter passed to ``synthesizeNTF`` which is used to compute the
        NTF zero location and the initial pole location for the optimizer

    Other Parameters
    ----------------
    show_progress : bool, optional
        provide extended output, default is False
    slsqp_xxx : various type, optional
        Parameters prefixed by ``slsqp_`` are passed to the ``fmin_slsqp``
        optimizer. Allowed options are:

        ``slsqp_maxiters``
            Maximum number of iterations (defaults to 100)
        ``slsqp_ftol``
            Precision goal for the value of f in the stopping criterion
            (defaults to 1e-6)
        ``slsqp_eps``
            Step size used for numerical approximation of the jacobian
            (defaults to 1.4901161193847656e-08)

        Do not use other options since they could break the minimizer in
        unexpected ways. Defaults can
        be set by changing the function ``default_options`` attribute.

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

    Notes
    -----
    This is actually an adaptation of the CLANS method. In the
    original paper about CLANS, NTF zeros are placed at z=1 (with the
    indication that other placements are possible). Here, the zeros are
    placed where ``synthesizeNTF`` places them.

    The computation is based on a nonlinear, nonlinearly constrained
    optimization. Since the optimizer used here is different from the
    optimizer used in other toolboxes implementing this funciton,
    the results may differ. The current optimizer is ``SLSQP``.

    The function internally calls ``synthesizeNTF``, and rises the same
    exceptions as ``synthesizeNTF``.
    """
    # Manage optional parameters
    opts = clans.default_options.copy()
    opts.update(options)
    slsqp_opts = {k[6:]: v for k, v in opts.iteritems()
                  if k.startswith('slsqp_')}
    if 'show_progress' in opts:
        slsqp_opts['disp'] = opts['show_progress']
    # Create the initial guess
    NTF = synthesizeNTF(order, osr, opt, 1+nq, 0)
    Hz = NTF[0]
    # x is the vector to optimize, accounting for the poles of the NTF
    # Rather than putting in x the poles of the NTF, the NTF is converted
    # into a CT form via a bilinear transform, then the frequency and damping
    # factor of each pole couple is recorded. This allows the constraint
    # that complex poles must be conjugated in couples to be automatically
    # managed. If the number of poles is odd, then there is necessarily a pole
    # that remains single, which should be treated differently, by storing
    # its frequency only.
    # Some notes...
    # 1. When practicing the bilinear transform to pass into the s-domain
    #    rather than using the poles of the NTF, one uses a version of such
    #    poles scaled by 1/rmax
    # 2. Rather than storing the poles frequency and damping factor, their
    #    square root is saved
    x = np.zeros(order)
    odd = order % 2
    poles = cplxpair(NTF[1])
    poles = poles[::-1]
    if odd:
        # Treat first pole (that is necessarily real)
        z = poles[0].real/rmax
        if np.abs(z) > 1:
            # Project poles outside rmax onto the circle with radius rmax
            # NB. Test has been fixed wrt DELSIG, where it always fails.
            z = z/np.abs(z)
        # Bilinear transform (modulus a mult. coeff)
        s = (z-1)/(z+1)
        x[0] = np.sqrt(-s)
    # All other poles are assumed to be complex.
    # N.B. This may be a wrong assumption, though. There may be more real
    # poles.
    for i in xrange(odd, order, 2):
        z = poles[i:i+2]/rmax
        if np.any(np.abs(z) > 1):
            # Project poles outside rmax onto the circle with radius rmax
            z = z/np.abs(z)
        # Bilinear transform (modulus a mult. coeff)
        s = (z-1)/(z+1)
        coeffs = np.poly(s).real
        wn = np.sqrt(coeffs[2])
        zeta = coeffs[1]/(2*wn)
        x[i] = np.sqrt(zeta)
        x[i+1] = np.sqrt(wn)

    # Run the optimizer
    x = minimize(_dsclansObj6a, x, args=(order, osr, nq, rmax, Hz),
                 method='SLSQP',
                 constraints={'type': 'ineq',
                              'fun': _dsclansObj6b,
                              'args': (order, osr, nq, rmax, Hz)},
                 options=slsqp_opts)['x']
    return dsclansNTF(x, order, rmax, Hz)

clans.default_options = {'show_progress': False,
                         'slsqp_maxiters': 100,
                         'slsqp_ftol': 1e-06,
                         'slsqp_eps': 1.4901161193847656e-08}


def _dsclansObj6a(x, order, osr, nq, rmax, Hz):
    # This is the objective function for clans
    # it returns the magnitude of H at the band edge.
    H = dsclansNTF(x, order, rmax, Hz)
    f = np.abs(evalTF(H, np.exp(1j*np.pi/osr)))
    return f


def _dsclansObj6b(x, order, osr, nq, rmax, Hz):
    # This is the constraint function for clans
    # it returns ||h||_1 - nq
    H = dsclansNTF(x, order, rmax, Hz)
    # NB. In the following use dB=60 rather than m=100, maybe
    g = np.sum(np.abs(impulse_response(H, m=100)))-1-nq
    # With our optimizer, this needs to be inverted
    return -g
