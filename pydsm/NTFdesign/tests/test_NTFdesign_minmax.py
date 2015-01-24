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

from __future__ import division, print_function

from numpy.testing import TestCase, run_module_suite
import numpy as np

from pydsm.NTFdesign import ntf_fir_minmax
from pydsm.relab import cplxpair

__all__ = ["TestSynthesizeNTFminmax"]


class TestSynthesizeNTFminmax(TestCase):
    def setUp(self):
        pass

    def test_LP8(self):
        z, p, k = ntf_fir_minmax(order=8, show_progress=False)
        e_k = 1
        e_z = [990.349427225477e-003 + 69.0500612157020e-003j,
               990.349427225477e-003 - 69.0500612157020e-003j,
               166.532844346146e-003 + 591.251073811726e-003j,
               166.532844346146e-003 - 591.251073811726e-003j,
               -259.915617496087e-003 + 503.342225950477e-003j,
               -259.915617496087e-003 - 503.342225950477e-003j,
               -512.031157651993e-003 + 194.699627385223e-003j,
               -512.031157651993e-003 - 194.699627385223e-003j]
        e_z = cplxpair(e_z)
        z = cplxpair(z)
        np.testing.assert_allclose(k, e_k, rtol=1e-6)
        np.testing.assert_allclose(z, e_z, rtol=2e-4)

    def test_BP8(self):
        z, p, k = ntf_fir_minmax(order=8, osr=32, f0=0.2, show_progress=False)
        e_k = 1
        e_z = [2.94348009789963e-01 + 9.14543800193135e-01j,
               2.94348009789963e-01 - 9.14543800193135e-01j,
               6.76745367518838e-01 + 0.00000000000000e+00j,
               2.46816733211163e-01 + 5.50000475735513e-01j,
               2.46816733211163e-01 - 5.50000475735513e-01j,
               -4.58884378359569e-01 + 4.10643263860101e-01j,
               -4.58884378359569e-01 - 4.10643263860101e-01j,
               -5.91022020183929e-01 + 0.00000000000000e+00j]
        e_z = cplxpair(e_z)
        z = cplxpair(z)
        np.testing.assert_allclose(k, e_k, rtol=1e-6)
        np.testing.assert_allclose(z, e_z, rtol=2e-4)


if __name__ == '__main__':
    run_module_suite()
