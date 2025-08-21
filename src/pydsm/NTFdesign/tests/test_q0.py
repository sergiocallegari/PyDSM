# -*- coding: utf-8 -*-

# Copyright (c) 2012–2024, Sergio Callegari
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

from __future__ import division, print_function

import numpy as np
from scipy import signal
from pydsm.ir import impulse_response
from pydsm.NTFdesign.legacy import q0_from_filter_ir
from pydsm.NTFdesign.weighting import q0_weighting
from pydsm.NTFdesign import quantization_noise_gain
import scipy.linalg as la

__all__ = ["TestQ0"]


class TestQ0:

    def setUp(self):
        pass

    def test_q0_butt_bp8(self):
        # Generate filter.
        # 8th order bandpass filter
        # Freq. passed to butterworth is normalized between 0 and 1
        # where 1 is the Nyquist frequency
        fsig = 1000.
        B = 400.
        OSR = 64
        fphi = B*OSR*2
        w0 = 2*fsig/fphi
        B0 = 2*B/fphi
        w1 = (np.sqrt(B0**2+4*w0**2)-B0)/2
        w2 = (np.sqrt(B0**2+4*w0**2)+B0)/2
        hz = signal.butter(4, [w1, w2], 'bandpass', output='zpk')
        # Order
        P = 20
        # Compute q0 in two ways
        ir = impulse_response(hz, db=80)

        q0_ir = q0_from_filter_ir(P, ir)
        q0_mr = q0_weighting(P, hz)
        np.testing.assert_allclose(q0_ir, q0_mr, atol=1E-7, rtol=1E-5)

    def test_q0_butt_lp3(self):
        # Generate filter.
        # 3rd order lowpass filter
        # Freq. passed to butterworth is normalized between 0 and 1
        # where 1 is the Nyquist frequency
        B = 400.
        OSR = 256
        fphi = B*OSR*2
        w0 = 2*B/fphi
        hz = signal.butter(3, w0, 'lowpass', output='zpk')
        # Order
        P = 12
        # Compute q0 in two ways
        ir = impulse_response(hz, db=80)

        q0_ir = q0_from_filter_ir(P, ir)
        q0_mr = q0_weighting(P, hz)
        np.testing.assert_allclose(q0_ir, q0_mr, atol=1E-7, rtol=1E-5)

    def test_q0_equiv(self):
        fir = np.asarray([1.00000000e+00,  -7.63387347e-01,  -4.02004111e-01,
                          -1.53885083e-01,  -2.76434316e-04,   7.91252937e-02,
                          1.03165832e-01,   8.83276154e-02,   4.92894661e-02])
        hz = (np.asarray([-1., -1., -1.]),
              np.asarray([0.99382681+0.01056265j,  0.98780284+0.j,
                          0.99382681-0.01056265j]),
              2.2820568419526305e-07)
        gain1 = quantization_noise_gain((fir, np.hstack((1, np.zeros(8)))), hz)
        q0 = q0_weighting(8, hz)
        Q = la.toeplitz(q0)
        fir_coeff = fir.reshape((1, 9))
        gain2 = fir_coeff.dot(Q).dot(fir_coeff.T)
        np.testing.assert_allclose(gain2, gain1)
