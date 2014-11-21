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

import unittest
import numpy as np
import scipy as sp
__import__("scipy.signal")
from pydsm.ir import impulse_response
from pydsm.delsig import evalTF
from pydsm.NTFdesign.filter_based import q0_from_filter

__all__=["TestQ0"]

class TestQ0(unittest.TestCase):

    def setUp(self):
        pass

    def test_q0_butt_bp8(self):
        # Generate filter.
        # 8th order bandpass filter
        # Freq. passed to butterworth is normalized between 0 and 1
        # where 1 is the Nyquist frequency
        fsig=1000.
        B=400.
        OSR=64
        fphi=B*OSR*2
        w0=2*fsig/fphi
        B0=2*B/fphi
        w1=(np.sqrt(B0**2+4*w0**2)-B0)/2
        w2=(np.sqrt(B0**2+4*w0**2)+B0)/2
        hz=sp.signal.butter(4, [w1,w2], 'bandpass', output='zpk')
        # Order
        P=20
        # Compute q0 in two ways
        ir=impulse_response(hz,db=80)
        def mr(f):
            return np.abs(evalTF(hz,np.exp(2j*np.pi*f)))
        q0_ir=q0_from_filter(P, ir, 'imp')
        q0_mr=q0_from_filter(P, mr, 'mag')
        np.testing.assert_allclose(q0_ir, q0_mr, atol=1E-7, rtol=1E-5)

    def test_q0_butt_lp3(self):
        # Generate filter.
        # 3rd order lowpass filter
        # Freq. passed to butterworth is normalized between 0 and 1
        # where 1 is the Nyquist frequency
        B=400.
        OSR=256
        fphi=B*OSR*2
        w0=2*B/fphi
        hz=sp.signal.butter(3, w0, 'lowpass', output='zpk')
        # Order
        P=12
        # Compute q0 in two ways
        ir=impulse_response(hz,db=80)
        def mr(f):
            return np.abs(evalTF(hz,np.exp(2j*np.pi*f)))
        q0_ir=q0_from_filter(P, ir, 'imp')
        q0_mr=q0_from_filter(P, mr, 'mag')
        np.testing.assert_allclose(q0_ir, q0_mr, atol=1E-7, rtol=1E-5)

if __name__ == '__main__':
    unittest.main()
