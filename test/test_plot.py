# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

import unittest

from pydsm.delsig import synthesizeNTF
from pydsm.delsig import plotPZ

__all__=["TestPlotPZ"]

class TestPlotPZ(unittest.TestCase):

    def setUp(self):
        pass
    
    def test_default(self):
        ntf = synthesizeNTF()
        plotPZ(ntf)

if __name__ == '__main__':
    unittest.main()
