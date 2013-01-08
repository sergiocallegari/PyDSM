# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

import unittest
import numpy as np
import csv
from pydsm.delsig import simulateDSM

__all__=["TestSimulateDSM"]

class TestSimulateDSM(unittest.TestCase):

    def setUp(self):
        pass

    def test_default(self):
        f = open('test/Data/test_simulateDSM_0.csv','r')
        csv_lines = csv.reader(f)
        d=np.array(csv_lines.next(),int)
        f.close()
        # Take H as in H = synthesizeNTF(5, 32, 1)
        H = (np.array([ 0.99604531+0.08884669j,  0.99604531-0.08884669j,
                    0.99860302+0.05283948j,  0.99860302-0.05283948j,
                    1.00000000+0.j ]),
             np.array([ 0.80655696+0.11982271j,  0.80655696-0.11982271j,
                    0.89807098+0.21981939j,  0.89807098-0.21981939j,
                    0.77776708+0.j ]),
             1)
        N = 8192
        f=85
        u = 0.5*np.sin(2.*np.pi*f/N*np.arange(N))
        v, d1, d2, d3 = simulateDSM(u, H)
        np.testing.assert_equal(v, d)

if __name__ == '__main__':
    unittest.main()
