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

import numpy as np
import time
import warnings
from pkg_resources import resource_stream
from numpy.testing import dec
from nose.plugins.skip import SkipTest

__all__ = ["Bench_simulateDSM"]


class Bench_simulateDSM():

    @classmethod
    def setup_class(cls):
        # Load expected results
        f = resource_stream('pydsm.delsig',
                            'benchmarks/Data/bench_simulateDSM_0.npz')
        cls.result = np.load(f)['arr_0']
        f.close()
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

    @classmethod
    def teardown_class(cls):
        pass

    def bench_simulateDSM_scipy_blas(self):
        """Benchmark function for the scipy blas version of simulateDSM"""
        from pydsm.delsig._simulateDSM_scipy_blas import (
            simulateDSM as simulateDSM_scipyblas)
        print("Benchmarking simulateDSM with scipy blas libraries")
        tic = time.clock()
        output, da1, da2, da3 = simulateDSM_scipyblas(self.u, self.H)
        timing = time.clock()-tic
        np.testing.assert_equal(output, self.result)
        print("Scipy Blas DSM simulator timing: %6.2f" % timing)

    def bench_simulateDSM_cblas_blas(self):
        """Benchmark function for the cblas version of simulateDSM"""
        try:
            from pydsm.delsig._simulateDSM_cblas import (
                simulateDSM as simulateDSM_cblas)
        except:
            raise SkipTest("Cblas libraries not available")
        print("Benchmarking simulateDSM with cblas libraries")
        tic = time.clock()
        output, da1, da2, da3 = simulateDSM_cblas(self.u, self.H)
        timing = time.clock()-tic
        np.testing.assert_equal(output, self.result)
        print("CBlas DSM simulator timing: %6.2f" % timing)

    @dec.slow
    def bench_simulateDSM_scipy(self):
        """Benchmark function for the scipy version of simulateDSM"""
        from pydsm.delsig._simulateDSM_scipy import (
            simulateDSM as simulateDSM_scipy)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            print("Benchmarking simulateDSM with scipy code")
            tic = time.clock()
            output, da1, da2, da3 = simulateDSM_scipy(self.u, self.H)
            timing = time.clock()-tic
        np.testing.assert_equal(output, self.result)
        print("Scipy DSM simulator timing: %6.2f" % timing)
