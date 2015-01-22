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

u"""
MinMax FIR NTF design (:mod:`pydsm.NTFdesign.minmax`)
=====================================================

This modules provides a minmax FIR design method for the NTF
of ΔΣ modulators

.. currentmodule:: pydsm.NTFdesign.minmax


Functions
---------

.. autosummary::
   :toctree: generated/

   ntf_fir_minmax    -- Design FIR NTF with minmax approach


Deprecated Functions
--------------------

.. autosummary::
   :toctree: generated/

   synthesize_ntf_minmax -- Alias of ntf_fir_minmax


Notes
-----
The design strategy implemented in this module was originally proposed by
M. Nagahara and Y. Yamamoto. See the documentation of :func:`ntf_fir_minmax`
for further details.
"""

from __future__ import division, print_function

import numpy as np
from scipy.signal import ss2zpk
import cvxpy_tinoco
from warnings import warn
from ..exceptions import PyDsmDeprecationWarning
from ..utilities import digested_options

__all__ = ['ntf_fir_minmax', 'synthesize_ntf_minmax']


def ntf_fir_minmax(order=32, osr=32, H_inf=1.5, f0=0, zf=False,
                   **options):
    u"""
    Synthesize FIR NTF for LP or BP ΔΣ modulator by min-max optimization.

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
        Provide extended output. Default is True and can be updated by
        changing the function ``default_options`` attribute.
    cvxopt_opts : dictionary, optional
        A dictionary of options for the ``cvxopt`` optimizer
        Allowed options include:

        ``maxiters``
            Maximum number of iterations (defaults to 100)
        ``abstol``
            Absolute accuracy (defaults to 1e-7)
        ``reltol``
            Relative accuracy (defaults to 1e-6)
        ``feastol``
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break ``cvxpy`` in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    Notes
    -----
    The design strategy implemented in this module is described in the paper
    [1]_.

    Bandpass modulator design is not yet supported.

    .. [1] M. Nagahara and Y. Yamamoto, *Frequency-Domain Min-Max Optimization
       of Noise-Shaping Delta-Sigma Modulators*, IEEE Trans. SP, vol. 60 n. 6
       June 2012.

    See also
    --------
    cvxopt : for the optimizer parameters.
    """
    # Manage optional parameters
    opts = digested_options(options, ntf_fir_minmax.default_options,
                            ['show_progress'], ['cvxopt_opts'])
    quiet = not opts['show_progress']

    # Maximum signal bandwidth in angular frequency
    Omega = 1./osr*np.pi
    # Centerband in angular frequency
    omega0 = 2*f0*np.pi

    # State space representation of NTF
    A = cvxpy_tinoco.matrix(
        np.vstack((np.hstack((np.zeros((order-1, 1)), np.eye(order-1))),
                   np.zeros((1, order)))))
    B = cvxpy_tinoco.matrix(np.vstack((np.zeros((order-1, 1)), [1.])))
    # C contains the NTF coefficients
    D = cvxpy_tinoco.matrix([[1]])

    # Set up the problem
    c = cvxpy_tinoco.variable(1, order)
    P = cvxpy_tinoco.variable(order, order, 'symmetric')
    Q = cvxpy_tinoco.variable(order, order, 'symmetric')
    g = cvxpy_tinoco.variable(1, 1)
    F = []
    if f0 == 0:
        # Lowpass modulator
        M1 = A.T*P*A+Q*A+A.T*Q-P-2*Q*np.cos(Omega)
        M2 = A.T*P*B + Q*B
        M3 = B.T*P*B-g
        M = cvxpy_tinoco.vstack(
            (cvxpy_tinoco.hstack((M1, M2, c.T)),
             cvxpy_tinoco.hstack((M2.T, M3, D)),
             cvxpy_tinoco.hstack((c, D, -1))))
        F += [cvxpy_tinoco.belongs(Q, cvxpy_tinoco.semidefinite_cone)]
        F += [cvxpy_tinoco.belongs(-M, cvxpy_tinoco.semidefinite_cone)]
        if zf:
            # Force a zero at DC
            F += [cvxpy_tinoco.equals(cvxpy_tinoco.sum(c), -1)]
    else:
        # Bandpass modulator
        M1r = (A.T*P*A + Q*A*np.cos(omega0) + A.T*Q*np.cos(omega0) -
               P - 2*Q*np.cos(Omega))
        M2r = A.T*P*B + Q*B*np.cos(omega0)
        M3r = B.T*P*B - g
        M1i = A.T*Q*np.sin(omega0) - Q*A*np.sin(omega0)
        M21i = -Q*B*np.sin(omega0)
        M22i = B.T*Q*np.sin(omega0)
        Mr = cvxpy_tinoco.vstack(
            (cvxpy_tinoco.hstack((M1r, M2r, c.T)),
             cvxpy_tinoco.hstack((M2r.T, M3r, D)),
             cvxpy_tinoco.hstack((c, D, -1))))
        Mi = cvxpy_tinoco.vstack(
            (cvxpy_tinoco.hstack((M1i, M21i, cvxpy_tinoco.zeros((order, 1)))),
             cvxpy_tinoco.hstack((M22i, 0, 0)),
             cvxpy_tinoco.hstack((cvxpy_tinoco.zeros((1, order)), 0, 0))))
        M = cvxpy_tinoco.vstack(
            (cvxpy_tinoco.hstack((Mr, Mi)),
             cvxpy_tinoco.hstack((-Mi, Mr))))
        F += [cvxpy_tinoco.belongs(Q, cvxpy_tinoco.semidefinite_cone)]
        F += [cvxpy_tinoco.belongs(-M, cvxpy_tinoco.semidefinite_cone)]
        if zf:
            # Force a zero at z=np.exp(1j*omega0)
            nn = np.arange(order).reshape((order, 1))
            vr = cvxpy_tinoco.matrix(np.cos(omega0*nn))
            vi = cvxpy_tinoco.matrix(np.sin(omega0*nn))
            vn = cvxpy_tinoco.matrix(
                [-np.cos(omega0*order), -np.sin(omega0*order)])
            F += [cvxpy_tinoco.equals(c*cvxpy_tinoco.hstack((vr, vi)), vn)]
    if H_inf < np.inf:
        # Enforce the Lee constraint
        R = cvxpy_tinoco.variable(order, order, 'symmetric')
        F += [cvxpy_tinoco.belongs(R, cvxpy_tinoco.semidefinite_cone)]
        MM = cvxpy_tinoco.vstack(
            (cvxpy_tinoco.hstack((A.T*R*A-R, A.T*R*B, c.T)),
             cvxpy_tinoco.hstack((B.T*R*A, -H_inf**2+B.T*R*B, D)),
             cvxpy_tinoco.hstack((c, D, -1))))
        F += [cvxpy_tinoco.belongs(-MM, cvxpy_tinoco.semidefinite_cone)]
    F += [cvxpy_tinoco.greater_equals(g, 0)]
    p = cvxpy_tinoco.program(cvxpy_tinoco.minimize(g), F)
    p.options.update(opts["cvxopt_opts"])
    p.solve(quiet)
    ntf = ss2zpk(A, B, np.asarray(c.value), D)
    return ntf

ntf_fir_minmax.default_options = {"cvxopt_opts": {'maxiters': 100,
                                                  'abstol': 1e-7,
                                                  'reltol': 1e-6,
                                                  'feastol': 1e-6},
                                  'show_progress': True}


# Following part is deprecated


def synthesize_ntf_minmax(order=32, osr=32, H_inf=1.5, f0=0, zf=False,
                          **options):
    """
    Alias of :func:`ntf_fir_minmax`

    .. deprecated:: 0.11.0
       Function is now available from the :mod:`NTFdesign` module with
       name :func:`ntf_fir_minmax`
    """
    warn("Function superseded by ntf_fir_minmax in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_fir_minmax(order, osr, H_inf, f0, zf, **options)
