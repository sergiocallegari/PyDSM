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

from __future__ import division, print_function

from numpy.testing import TestCase, run_module_suite
import numpy as np
from scipy import signal
from pydsm.ir import impulse_response
from pydsm.NTFdesign import (ntf_fir_weighting, quantization_noise_gain,
                             mult_weightings)
from pydsm.NTFdesign.weighting import ntf_fir_from_q0
from pydsm.NTFdesign.legacy import (quantization_noise_gain_by_conv,
                                    q0_from_filter_ir)
from nose.plugins.skip import SkipTest

__all__ = ["TestNTF_Filter"]


class TestNTF_Filter(TestCase):

    def setUp(self):
        self.z_e = np.asarray([0.98979462+0.12667657j,
                               0.98979462-0.12667657j,
                               0.72084151+0.0j,
                               0.35347507+0.64857031j,
                               0.35347507-0.64857031j,
                               -0.02875404+0.70480695j,
                               -0.02875404-0.70480695j,
                               -0.36294495+0.58281858j,
                               -0.36294495-0.58281858j,
                               -0.67350105+0.0j,
                               -0.59201143+0.32765994j,
                               -0.59201143-0.32765994j])
        self.z_e = np.sort(self.z_e)
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
        self.hz = signal.butter(4, [w1, w2], 'bandpass', output='zpk')
        self.order = self.z_e.size

    def test_ntf_butt_bp8_vs_legacy(self):
        # Compute q0 in two ways
        ir = impulse_response(self.hz, db=80)
        ntf1 = ntf_fir_weighting(self.order, self.hz, show_progress=False)
        q0 = q0_from_filter_ir(self.order, ir)
        ntf2 = ntf_fir_from_q0(q0, show_progress=False)
        mf1 = quantization_noise_gain(ntf1, self.hz)
        mf2 = quantization_noise_gain(ntf2, self.hz)
        np.testing.assert_allclose(mf2, mf1, rtol=1e-7)
        mf3 = quantization_noise_gain_by_conv(ntf1, self.hz)
        np.testing.assert_allclose(mf3, mf1, rtol=5e-6)

    def test_ntf_butt_bp8_cvxpy_old(self):
        try:
            import cvxpy_tinoco     # analysis:ignore
        except:
            raise SkipTest("Modeler 'cvxpy_old' not installed")
        ntf1 = ntf_fir_weighting(self.order, self.hz, modeler='cvxpy_old',
                                 show_progress=False)
        np.testing.assert_allclose(np.sort(ntf1[0]), self.z_e, rtol=1e-5)

    def test_ntf_butt_bp8_cvxpy_cvxopt(self):
        try:
            import cvxpy     # analysis:ignore
        except:
            raise SkipTest("Modeler 'cvxpy' not installed")
        ntf1 = ntf_fir_weighting(self.order, self.hz, modeler='cvxpy',
                                 show_progress=False)
        np.testing.assert_allclose(np.sort(ntf1[0]), self.z_e, rtol=1e-7)

    def test_ntf_butt_bp8_cvxpy_scs(self):
        try:
            import cvxpy     # analysis:ignore
        except:
            raise SkipTest("Modeler 'cvxpy' not installed")
        ntf1 = ntf_fir_weighting(self.order, self.hz, modeler='cvxpy',
                                 cvxpy_opts={'solver': 'scs'},
                                 show_progress=False)
        np.testing.assert_allclose(np.sort(ntf1[0]), self.z_e, rtol=1e-4)

    def test_ntf_butt_bp8_picos(self):
        try:
            import picos     # analysis:ignore
        except:
            raise SkipTest("Modeler 'picos' not installed")
        ntf1 = ntf_fir_weighting(self.order, self.hz, modeler='picos',
                                 show_progress=False)
        np.testing.assert_allclose(np.sort(ntf1[0]), self.z_e, rtol=1e-7)


class Test_MultWeightings(TestCase):

    def setUp(self):
        pass

    def test_mult_weightings(self):
        f = mult_weightings(([], [], 1), ([], [], 2))
        np.testing.assert_equal(f(0), 4)

if __name__ == '__main__':
    run_module_suite()
