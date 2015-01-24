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

    def test_MB16(self):
        z, p, k = ntf_fir_minmax(order=16, osr=64, f0=[0.1, 0.2],
                                 show_progress=False)
        e_k = 1
        e_z = [7.987069488919752e-01 + 5.788064442124057e-01j,
               7.987069488919752e-01 - 5.788064442124057e-01j,
               3.043919962928907e-01 + 9.396829116486076e-01j,
               3.043919962928907e-01 - 9.396829116486076e-01j,
               7.613233234492194e-01,
               4.324139559044125e-01 + 3.504979544039063e-01j,
               4.324139559044125e-01 - 3.504979544039063e-01j,
               -6.539374591580128e-01,
               -6.069758142797348e-01 + 2.496966124458908e-01j,
               -6.069758142797348e-01 - 2.496966124458908e-01j,
               -4.679294834745639e-01 + 4.717700769351812e-01j,
               -4.679294834745639e-01 - 4.717700769351812e-01j,
               -2.310284089771453e-01 + 6.478900732051520e-01j,
               -2.310284089771453e-01 - 6.478900732051520e-01j,
               -5.901531897720751e-02 + 4.677138165783495e-01j,
               -5.901531897720751e-02 - 4.677138165783495e-01j]
        e_z = cplxpair(e_z)
        z = cplxpair(z)
        np.testing.assert_allclose(k, e_k, rtol=1e-6)
        np.testing.assert_allclose(z, e_z, rtol=1e-1)

if __name__ == '__main__':
    run_module_suite()
