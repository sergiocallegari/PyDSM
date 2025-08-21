# -*- coding: utf-8 -*-

# Copyright (c) 2012-2024, Sergio Callegari
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
from pydsm.delsig import evalTF
from pydsm.relab import db

__all__ = ["TestIR"]


class TestIR:

    def setUp(self):
        pass

    def test_butt_bp8_ir(self):
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
        ir = impulse_response(hz, db=80)
        ff = np.logspace(np.log10(w0/4), np.log10(w0*4), 100)
        vv1 = evalTF(hz, np.exp(2j*np.pi*ff))
        vv2 = evalTF((ir, [1]), np.exp(2j*np.pi*ff))
        np.testing.assert_allclose(db(vv1/vv2), np.zeros_like(vv1),
                                   rtol=0, atol=3)

    def test_fir_ir(self):
        fir = np.asarray([1.0, 0.5, 0.25, 0.125])
        zpk = (np.roots(fir), np.zeros(4), 1)
        ir = impulse_response(zpk, db=80)
        np.testing.assert_allclose(fir, ir)
