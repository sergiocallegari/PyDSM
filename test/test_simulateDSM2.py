# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

import unittest
import numpy as np
skip_test=False
try:
    from pydsm.delsig._simulateDSM_cblas import simulateDSM
except:
    skip_test=True
from pydsm.delsig._simulateDSM_scipy_blas import simulateDSM as simulateDSM2
import time

__all__=["TestSimulateDSM2"]

class TestSimulateDSM2(unittest.TestCase):

    def setUp(self):
        pass

    def test_default(self):
        if skip_test:
            print "Skipping on this platform"
            return
        # Take H as in H = synthesizeNTF(5, 32, 1)
        H = (np.array([ 0.99604531+0.08884669j,  0.99604531-0.08884669j,
                    0.99860302+0.05283948j,  0.99860302-0.05283948j,
                    1.00000000+0.j ]),
             np.array([ 0.80655696+0.11982271j,  0.80655696-0.11982271j,
                    0.89807098+0.21981939j,  0.89807098-0.21981939j,
                    0.77776708+0.j ]),
             1)
        N = 819200
        f=85
        u = 0.5*np.sin(2.*np.pi*f/N*np.arange(N))
        tic=time.clock()
        va, da1, da2, da3 = simulateDSM2(u, H)
        tac=time.clock()
        vb, db1, db2, db3 = simulateDSM(u, H)
        toc=time.clock()
        np.testing.assert_equal(va, vb)
        print "[scipy blas: ", tac-tic, "cblas: ", toc-tac, "]"

if __name__ == '__main__':
    unittest.main()
