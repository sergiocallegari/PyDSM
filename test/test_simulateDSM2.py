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
