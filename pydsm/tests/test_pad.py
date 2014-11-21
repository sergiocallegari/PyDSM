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
