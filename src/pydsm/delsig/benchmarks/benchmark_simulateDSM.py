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
import warnings
import importlib_resources
import pytest

try:
    import pytest_benchmark
    BENCHMARK_AVAILABLE = True
except ImportError:
    BENCHMARK_AVAILABLE = False

@pytest.mark.skipif(not BENCHMARK_AVAILABLE,
                    reason="pytest-benchmark is not installed")
@pytest.mark.benchmark(group="simulator")
class Benchmark_simulateDSM:

    @classmethod
    def setup_class(cls):
        # Load expected results
        with (importlib_resources.files('pydsm.delsig')
              .joinpath('tests/Data/test_simulateDSM_0.npz')
              .open('rb')) as f:
            cls.result = np.load(f)['arr_0']
        cls.N = cls.result.size
        # Take H as in H = synthesizeNTF(5, 32, 1)
        cls.H = (np.array([0.99604531+0.08884669j,  0.99604531-0.08884669j,
                           0.99860302+0.05283948j,  0.99860302-0.05283948j,
                           1.00000000+0.j]),
                 np.array([0.80655696+0.11982271j,  0.80655696-0.11982271j,
                           0.89807098+0.21981939j,  0.89807098-0.21981939j,
                           0.77776708+0.j]),
                 1)
        # Test frequency
        cls.f = 85.
        # Test input signal
        cls.u = 0.5*np.sin(2.*np.pi*cls.f/cls.N*np.arange(cls.N))

    def benchmark_simulateDSM_scipy_blas(self, benchmark):
        """Benchmark function for the scipy blas version of simulateDSM"""
        from pydsm.delsig._simulateDSM_scipy_blas import (
            simulateDSM as simulateDSM_scipyblas)
        output, da1, da2, da3 = benchmark(simulateDSM_scipyblas,
                                          self.u, self.H)
        np.testing.assert_equal(output, self.result)

    def benchmark_simulateDSM_cblas_blas(self, benchmark):
        """Benchmark function for the cblas version of simulateDSM"""
        try:
            from pydsm.delsig._simulateDSM_cblas import (
                simulateDSM as simulateDSM_cblas)
        except:
            pytest.skip("Cblas libraries not available")
        output, da1, da2, da3 = benchmark(simulateDSM_cblas,
                                          self.u, self.H)
        np.testing.assert_equal(output, self.result)

    @pytest.mark.slow
    def benchmark_simulateDSM_scipy(self, benchmark):
        """Benchmark function for the scipy version of simulateDSM"""
        from pydsm.delsig._simulateDSM_scipy import (
            simulateDSM as simulateDSM_scipy)
        from pydsm.exceptions import PyDsmSlowPathWarning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",
                category=PyDsmSlowPathWarning)
            output, da1, da2, da3 = benchmark(simulateDSM_scipy,
                                              self.u, self.H)
        np.testing.assert_equal(output, self.result)
