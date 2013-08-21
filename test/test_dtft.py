# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

# -*- coding: utf-8 -*-

import unittest
import numpy as np
from pydsm.ft import dtft, dtft_hermitian, idtft, idtft_hermitian

__all__=["TestDTFT"]

class TestDTFT(unittest.TestCase):

    def setUp(self):
        pass

    def test_dtft_re10(self):
        xx=np.arange(1,11)
        Ff=dtft(xx)
        xx2=idtft(Ff,np.arange(10))
        np.testing.assert_allclose(xx,xx2, atol=1E-12)

    def test_dtft_re10s1(self):
        xx=np.arange(1,11)
        Ff=dtft_hermitian(xx)
        xx2=idtft(Ff,np.arange(10))
        np.testing.assert_allclose(xx,xx2, atol=1E-12)

    def test_dtft_re10s2(self):
        xx=np.arange(1,11)
        Ff=dtft_hermitian(xx)
        xx2=idtft_hermitian(Ff,np.arange(10))
        np.testing.assert_allclose(xx,xx2, atol=1E-12)

if __name__ == '__main__':
    unittest.main()
