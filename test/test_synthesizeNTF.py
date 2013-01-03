# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

import unittest
import numpy as np

from pydsm.delsig import synthesizeNTF
from pydsm.utilities import cplxpair

__all__=["TestSynthesizeNTF"]

class TestSynthesizeNTF(unittest.TestCase):

    def setUp(self):
        pass

    def test_default(self):
        z, p, k = synthesizeNTF()
        e_k = 1
        e_z = np.ones(3)
        e_p = [0.7654 - 0.2793j, 0.7654 + 0.2793j, 0.6694]
        e_p = cplxpair(e_p)
        self.assertEqual(k, e_k)
        np.testing.assert_equal(z, e_z)
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
        self.assertEqual(k, e_k)
        np.testing.assert_almost_equal(z, e_z, 4)
        np.testing.assert_almost_equal(p, e_p, 4)        

    def test_opt1(self):
        z, p, k = synthesizeNTF(opt=1)
        e_k = 1
        e_z = [1.0000, 0.9993 - 0.0380j, 0.9993 + 0.0380j]
        e_z = cplxpair(e_z)                    
        e_p = [0.6692, 0.7652 - 0.2795j, 0.7652 + 0.2795j]
        e_p = cplxpair(e_p)
        self.assertEqual(k, e_k)
        np.testing.assert_almost_equal(z, e_z, 4)
        np.testing.assert_almost_equal(p, e_p, 4)
        
    def test_opt3(self):
        z, p, k = synthesizeNTF(opt=3)
        e_k = 1
        e_z = [1.0000, 0.9993 - 0.0382j, 0.9993 + 0.0382j]
        e_z = cplxpair(e_z)                    
        e_p = [0.6692, 0.7652 - 0.2795j, 0.7652 + 0.2795j]
        e_p = cplxpair(e_p)
        self.assertEqual(k, e_k)
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
        self.assertEqual(k, e_k)
        np.testing.assert_almost_equal(z, e_z, 4)
        np.testing.assert_almost_equal(p, e_p, 4)

if __name__ == '__main__':
    unittest.main()
