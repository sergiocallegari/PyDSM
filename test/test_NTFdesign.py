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

from pydsm.NTFdesign.delsig import (synthesizeNTF, synthesizeChebyshevNTF,
                                    clans)
from pydsm.utilities import cplxpair

__all__=["TestSynthesizeNTF", "TestSynthesizeChebyshevNTF",
         "TestClans"]

class TestSynthesizeNTF(unittest.TestCase):
    def setUp(self):
        pass

    def test_default(self):
        z, p, k = synthesizeNTF()
        e_k = 1
        e_z = np.ones(3)
        e_p = [0.7654 - 0.2793j, 0.7654 + 0.2793j, 0.6694]
        e_p = cplxpair(e_p)
        np.testing.assert_almost_equal(k, e_k, 6)
        np.testing.assert_almost_equal(z, e_z, 6)
        np.testing.assert_almost_equal(p, e_p, 4)

    def test_BP6(self):
        z, p, k = synthesizeNTF(order=6, f0=0.3)
        e_k = 1
        e_z = np.concatenate(((-0.3090 + 0.9511j)*np.ones(3),
                            (-0.3090 - 0.9511j)*np.ones(3)))
        e_z = cplxpair(e_z)
        e_p = [-0.4238 - 0.7906j, -0.4238 + 0.7906j, -0.2480 - 0.7802j,
               -0.2480 + 0.7802j, -0.1290 - 0.8985j, -0.1290 + 0.8985j]
        e_p = cplxpair(e_p)
        np.testing.assert_almost_equal(k, e_k, 6)
        np.testing.assert_almost_equal(z, e_z, 4)
        np.testing.assert_almost_equal(p, e_p, 4)

    def test_opt1(self):
        z, p, k = synthesizeNTF(opt=1)
        e_k = 1
        e_z = [1.0000, 0.9993 - 0.0380j, 0.9993 + 0.0380j]
        e_z = cplxpair(e_z)
        e_p = [0.6692, 0.7652 - 0.2795j, 0.7652 + 0.2795j]
        e_p = cplxpair(e_p)
        np.testing.assert_almost_equal(k, e_k, 6)
        np.testing.assert_almost_equal(z, e_z, 4)
        np.testing.assert_almost_equal(p, e_p, 4)

    def test_opt3(self):
        z, p, k = synthesizeNTF(opt=3)
        e_k = 1
        e_z = [1.0000, 0.9993 - 0.0382j, 0.9993 + 0.0382j]
        e_z = cplxpair(e_z)
        e_p = [0.6692, 0.7652 - 0.2795j, 0.7652 + 0.2795j]
        e_p = cplxpair(e_p)
        np.testing.assert_almost_equal(k, e_k, 6)
        np.testing.assert_almost_equal(z, e_z, 4)
        np.testing.assert_almost_equal(p, e_p, 4)

    def test_given_zeroz(self):
        g_z = [1.0000, 0.9993 - 0.0380j, 0.9993 + 0.0380j]
        z, p, k = synthesizeNTF(opt=g_z)
        e_k = 1
        e_z = [1.0000, 0.9993 - 0.0380j, 0.9993 + 0.0380j]
        e_z = cplxpair(e_z)
        e_p = [0.6692, 0.7652 - 0.2795j, 0.7652 + 0.2795j]
        e_p = cplxpair(e_p)
        np.testing.assert_almost_equal(k, e_k, 6)
        np.testing.assert_almost_equal(z, e_z, 4)
        np.testing.assert_almost_equal(p, e_p, 4)

class TestSynthesizeChebyshevNTF(unittest.TestCase):

    def setUp(self):
        pass

    def test_8_4_3(self):
        order = 8
        osr = 4
        h_inf = 3
        z, p, k = synthesizeChebyshevNTF(order, osr, 1, h_inf)
        e_k = 1
        e_z = [716.675171843707e-003 + 697.407125044473e-003j,
               716.675171843707e-003 - 697.407125044473e-003j,
               987.024523733349e-003 + 160.569578528934e-003j,
               987.024523733349e-003 - 160.569578528934e-003j,
               787.924470008031e-003 + 615.771897347195e-003j,
               787.924470008031e-003 - 615.771897347195e-003j,
               899.412094416429e-003 + 437.101686587289e-003j,
               899.412094416429e-003 - 437.101686587289e-003j]
        e_p = [626.054825241768e-003 + 122.563529859845e-003j,
               626.054825241768e-003 - 122.563529859845e-003j,
               603.812988962951e-003 + 354.087262468194e-003j,
               603.812988962951e-003 - 354.087262468194e-003j,
               583.053683660931e-003 + 552.205015613472e-003j,
               583.053683660931e-003 - 552.205015613472e-003j,
               599.974045233308e-003 + 710.003702019636e-003j,
               599.974045233308e-003 - 710.003702019636e-003j]
        e_z = cplxpair(e_z)
        e_p = cplxpair(e_p)
        np.testing.assert_almost_equal(k, e_k, 4)
        np.testing.assert_almost_equal(z, e_z, 4)
        np.testing.assert_almost_equal(p, e_p, 4)


class TestClans(unittest.TestCase):
    def setUp(self):
        pass

    def test_5_32_5_095_1(self):
        z, p, k = clans(5,32,5,0.95,1)
        e_k = 1
        e_z = [1.00000000000000e+000,
               998.603018798634e-003 + 52.8394818885978e-003j,
               998.603018798634e-003 - 52.8394818885978e-003j,
               996.045312843931e-003 + 88.8466924631155e-003j,
               996.045312843931e-003 - 88.8466924631155e-003j]
        e_p = [418.352260015271e-003,
               489.222137893601e-003 + 170.971685170035e-003j,
               489.222137893601e-003 - 170.971685170035e-003j,
               652.448968088646e-003 + 381.722332008184e-003j,
               652.448968088646e-003 - 381.722332008184e-003j]
        e_z = cplxpair(e_z)
        e_p = cplxpair(e_p)
        z = cplxpair(z)
        p = cplxpair(p)
        np.testing.assert_almost_equal(k, e_k, 4)
        np.testing.assert_almost_equal(z, e_z, 4)
        np.testing.assert_almost_equal(p, e_p, 4)

if __name__ == '__main__':
    unittest.main()
