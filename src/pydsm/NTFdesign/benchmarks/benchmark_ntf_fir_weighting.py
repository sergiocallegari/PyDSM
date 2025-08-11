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
# along with PyDSM.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division, print_function

import numpy as np
from scipy import signal
from pydsm.NTFdesign import ntf_fir_weighting, quantization_noise_gain
from pydsm.delsig import evalTF
import pytest

try:
    import pytest_benchmark
    BENCHMARK_AVAILABLE = True
except ImportError:
    BENCHMARK_AVAILABLE = False

@pytest.mark.skipif(not BENCHMARK_AVAILABLE,
                    reason="pytest-benchmark is not installed")
@pytest.mark.benchmark(group="ntf_design-fir-weighting")
class Benchmark_ntf_fir_weighting(object):

    @classmethod
    def setup_class(cls):
        # signal and filter specification
        fsig = 1000.
        B = 400.
        OSR = 64
        fphi = B*OSR*2
        # Lee constraint
        cls.H_inf = 1.5
        # FIR Order
        cls.order = 25
        w0 = 2*fsig/fphi
        B0 = 2*B/fphi
        w1 = (np.sqrt(B0**2+4*w0**2)-B0)/2
        w2 = (np.sqrt(B0**2+4*w0**2)+B0)/2
        cls.hz = signal.butter(4, [w1, w2], 'bandpass', output='zpk')
        cls.ff = np.linspace(0, 0.5, 1024)

    @classmethod
    def teardown_class(cls):
        pass

    def benchmark_ntf_fir_weighting_cvxpy_old(self, benchmark):
        try:
            from pydsm import cvxpy_tdr     # analysis:ignore
        except:
            pytest.skip("Modeler 'cvxpy_old' not installed")
        ntf = benchmark(ntf_fir_weighting,
                        self.order, self.hz, self.H_inf,
                        modeler='cvxpy_old', show_progress=False)
        mf = quantization_noise_gain(ntf, self.hz)
        vv = np.abs(evalTF(ntf, np.exp(2j*np.pi*self.ff)))
        pe = np.max(vv)-self.H_inf
        benchmark.extra_info['Constraints delta'] = pe
        benchmark.extra_info['Quantization noise'] = mf

    def benchmark_ntf_fir_weighting_cvxpy_cvxopt(self, benchmark):
        try:
            import cvxpy     # analysis:ignore
        except:
            pytest.skip("Modeler 'cvxpy' not installed")
        ntf = benchmark(ntf_fir_weighting,
                        self.order, self.hz, self.H_inf,
                        modeler='cvxpy', show_progress=False)
        mf = quantization_noise_gain(ntf, self.hz)
        vv = np.abs(evalTF(ntf, np.exp(2j*np.pi*self.ff)))
        pe = np.max(vv)-self.H_inf
        benchmark.extra_info['Constraints delta'] = pe
        benchmark.extra_info['Quantization noise'] = mf

    def benchmark_ntf_fir_weighting_cvxpy_cvxopt_kkt(self, benchmark):
        try:
            import cvxpy     # analysis:ignore
        except:
            pytest.skip("Modeler 'cvxpy' not installed")
        ntf = benchmark(ntf_fir_weighting,
                        self.order, self.hz, self.H_inf,
                        modeler='cvxpy', show_progress=False,
                        cvxpy_opts={'override_kktsolver': False})
        mf = quantization_noise_gain(ntf, self.hz)
        vv = np.abs(evalTF(ntf, np.exp(2j*np.pi*self.ff)))
        pe = np.max(vv)-self.H_inf
        benchmark.extra_info['Constraints delta'] = pe
        benchmark.extra_info['Quantization noise'] = mf

    def benchmark_ntf_fir_weighting_cvxpy_scs(self, benchmark):
        try:
            import cvxpy     # analysis:ignore
        except:
            pytest.skip("Modeler 'cvxpy' not installed")
        ntf = benchmark(ntf_fir_weighting,
                        self.order, self.hz, self.H_inf,
                        modeler='cvxpy', show_progress=False,
                        cvxpy_opts={'solver': 'scs'})
        mf = quantization_noise_gain(ntf, self.hz)
        vv = np.abs(evalTF(ntf, np.exp(2j*np.pi*self.ff)))
        pe = np.max(vv)-self.H_inf
        benchmark.extra_info['Constraints delta'] = pe
        benchmark.extra_info['Quantization noise'] = mf

    def benchmark_ntf_fir_weighting_picos(self, benchmark):
        try:
            import picos     # analysis:ignore
        except:
            pytest.skip("Modeler 'picos' not installed")
        ntf = benchmark(ntf_fir_weighting,
                        self.order, self.hz, self.H_inf,
                        modeler='picos', show_progress=False)
        mf = quantization_noise_gain(ntf, self.hz)
        vv = np.abs(evalTF(ntf, np.exp(2j*np.pi*self.ff)))
        pe = np.max(vv)-self.H_inf
        benchmark.extra_info['Constraints delta'] = pe
        benchmark.extra_info['Quantization noise'] = mf
