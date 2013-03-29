# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

import unittest
import numpy as np

from pydsm.NTFdesign.delsig import synthesizeNTF, synthesizeChebyshevNTF
from pydsm.utilities import cplxpair

__all__=["TestSynthesizeNTF", "TestSynthesizeChebyshevNTF"]

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

if __name__ == '__main__':
    unittest.main()
