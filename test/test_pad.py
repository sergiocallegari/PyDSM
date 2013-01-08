# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

import unittest
import numpy as np

from pydsm.delsig import padl, padr, padt, padb

__all__=["TestSynthesizePAD"]

class TestSynthesizePAD(unittest.TestCase):

    def setUp(self):
        pass

    def test_padl_vector(self):
        v = np.asarray([1, 2, 3])
        ve = np.asarray([4, 4, 1, 2, 3])
        np.testing.assert_equal(padl(v, 5, 4), ve)

    def test_padr_vector(self):
        v = np.asarray([1, 2, 3])
        ve = np.asarray([1, 2, 3, 0, 0])
        np.testing.assert_equal(padr(v, 5), ve)

    def test_padl_array2d(self):
        v = np.asarray([[1, 2, 3], [4, 5, 6]])
        ve = np.asarray([[4, 4, 1, 2, 3],[4, 4, 4, 5, 6]])
        np.testing.assert_equal(padl(v, 5, 4), ve)

    def test_padr_array2d(self):
        v = np.asarray([[1, 2, 3], [4, 5, 6]])
        ve = np.asarray([[1, 2, 3, 0, 0],[4, 5, 6, 0, 0]])
        np.testing.assert_equal(padr(v, 5), ve)

    def test_padt_array2d(self):
        v = np.asarray([[1, 2, 3], [4, 5, 6]])
        ve = np.asarray([[4, 4, 4], [4, 4, 4], [1, 2, 3], [4, 5, 6]])
        np.testing.assert_equal(padt(v, 4, 4), ve)

    def test_padb_array2d(self):
        v = np.asarray([[1, 2, 3], [4, 5, 6]])
        ve = np.asarray([[1, 2, 3], [4, 5, 6], [0, 0, 0], [0, 0, 0]])
        np.testing.assert_equal(padb(v, 4), ve)

if __name__ == '__main__':
    unittest.main()
