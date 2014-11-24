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

# This file includes code ported from the "Frequency-Domain Min-Max
# Optimization for Delta-Sigma Modulators" Matlab toolbox by Masaaki
# Nagahara (see http://www.mathworks.com/matlabcentral/fileexchange/36187)
# covered by the following copyright and permission notice
#
# Copyright (c) 2012, Masaaki Nagahara
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the distribution
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""
Synthesize a FIR NTF using the min-max strategy by M. Nagahara et al
====================================================================

The design strategy implemented in this module is described in the paper

    M. Nagahara and Y. Yamamoto, *Frequency-Domain Min-Max Optimization of
    Noise-Shaping Delta-Sigma Modulators*, IEEE Trans. SP, 2012.
"""

from __future__ import division, print_function

import numpy as np
from scipy.signal import ss2zpk
import cvxpy_tinoco as cvxpy

__all__ = ['synthesize_ntf_minmax']


def synthesize_ntf_minmax(order=32, osr=32, H_inf=1.5, f0=0, zf=False,
                          **options):
    u"""Synthesize FIR NTF for LP or BP ΔΣ modulator by min-max optimization.

    The design strategy implemented in this module is described in the paper

        M. Nagahara and Y. Yamamoto, *Frequency-Domain Min-Max Optimization of
        Noise-Shaping Delta-Sigma Modulators*, IEEE Trans. SP, vol. 60 n. 6
        June 2012.

    The strategy aims at minimizing the peak value of the NTF in the signal
    band, while respecting the Lee criterion.

    Parameters
    ----------
    order : int, optional
        Order of the delta-sigma modulator. Defaults to 32.
    osr : real, optional
        The oversampling ratio, based on the input signal bandwidth
        Defaults to 32.
    H_inf : real, optional
        Max allowed peak value of the NTF. Used to enforce the Lee criterion.
        Defaults to 1.5.
    f0 : real, optional
        The normalized center frequency of the modulator. Value in [0,1/2].
        1/2 indicates half the sample frequency. Defaults to 0, indicating an
        LP modulator.
    zf : bool, optional
        Flag controlling the pre-assignement of NTF zeros. If ``False``, the
        design is practiced without any zero pre-assignment. If ``True``, a
        zero is pre-assigned at the modulator center-band. Defaults to False.

    Returns
    -------
    ntf : tuple
        noise transfer function in zpk form.

    Other parameters
    ----------------
    show_progress : bool, optional
        provide extended output, default is True
    cvxpy_xxx : various type, optional
        Parameters prefixed by ``cvxpy_`` are passed to the ``cvxpy``
        optimizer. Allowed options are:

        ``cvxpy_maxiters``
            Maximum number of iterations (defaults to 100)
        ``cvxpy_abstol``
            Absolute accuracy (defaults to 1e-7)
        ``cvxpy_reltol``
            Relative accuracy (defaults to 1e-6)
        ``cvxpy_feastol``
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break ``cvxpy`` in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    See Also
    --------

    Check the documentation of ``cvxopt`` for further information.

    Notes
    -----
    Bandpass modulator design is not yet supported.
    """
    # Manage optional parameters
    opts = synthesize_ntf_minmax.default_options.copy()
    opts.update(options)
    cvxpy_opts = {k[6:]: v for k, v in opts.iteritems()
                  if k.startswith('cvxpy_')}
    if 'show_progress' in opts:
        quiet = not opts['show_progress']
    else:
        quiet = False

    # Maximum signal frequency
    Omega = 1./osr*np.pi

    # State space representation of NTF
    A = cvxpy.matrix(np.vstack((np.hstack((np.zeros((order-1, 1)),
                                           np.eye(order-1))),
                                np.zeros((1, order)))))
    B = cvxpy.matrix(np.vstack((np.zeros((order-1, 1)), [1.])))
    # C contains the NTF coefficients
    D = cvxpy.matrix([[1]])

    # Set up the problem
    c = cvxpy.variable(1, order)
    P = cvxpy.variable(order, order, 'symmetric')
    Q = cvxpy.variable(order, order, 'symmetric')
    g = cvxpy.variable(1, 1)
    F = []
    if f0 == 0:
        # Lowpass modulator
        M1 = A.T*P*A+Q*A+A.T*Q-P-2*Q*np.cos(Omega)
        M2 = A.T*P*B + Q*B
        M3 = B.T*P*B-g
        M = cvxpy.vstack((cvxpy.hstack((M1, M2, c.T)),
                          cvxpy.hstack((M2.T, M3, D)),
                          cvxpy.hstack((c, D, -D))))
        F += [cvxpy.belongs(Q, cvxpy.semidefinite_cone)]
        F += [cvxpy.belongs(-M, cvxpy.semidefinite_cone)]
        if zf:
            # Force a zero at DC
            F += [cvxpy.equals(cvxpy.sum(c), -1)]
        if H_inf < np.inf:
            # Enforce the Lee constraint
            R = cvxpy.variable(order, order, 'symmetric')
            F += [cvxpy.belongs(R, cvxpy.semidefinite_cone)]
            MM = cvxpy.vstack((cvxpy.hstack((A.T*R*A-R, A.T*R*B, c.T)),
                               cvxpy.hstack((B.T*R*A, -H_inf**2+B.T*R*B, 1)),
                               cvxpy.hstack((c, D, -D))))
            F += [cvxpy.belongs(-MM, cvxpy.semidefinite_cone)]
    else:
        # Bandpass modulator
        # not implemented yet
        return None
    F += [cvxpy.greater_equals(g, 0)]
    p = cvxpy.program(cvxpy.minimize(g), F)
    p.options.update(cvxpy_opts)
    p.solve(quiet)
    ntf = ss2zpk(A, B, np.asarray(c.value), D)
    return ntf

synthesize_ntf_minmax.default_options = {}
