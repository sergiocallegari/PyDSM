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

from __future__ import division, print_function

import numpy as np
from numpy.testing import TestCase, run_module_suite
from pydsm.delsig import clans
from pydsm.utilities import cplxpair

__all__ = ["TestClans"]


class TestClans(TestCase):

    def setUp(self):
        pass

    def test_clans(self):
        """Test function for clans"""
        e_z = [1,
               0.99604531 + 0.08884669j,
               0.99604531 - 0.08884669j,
               0.99860302 + 0.05283948j,
               0.99860302 - 0.05283948j]
        e_p = [0.41835234,
               0.48922229 + 0.1709716j,
               0.48922229 - 0.1709716j,
               0.65244884 + 0.38172239j,
               0.65244884 - 0.38172239j]
        e_k = 1
        e_z = cplxpair(e_z)
        e_p = cplxpair(e_p)
        z, p, k = clans(5, 32, 5, .95, 1)
        z = cplxpair(z)
        p = cplxpair(p)
        np.testing.assert_almost_equal(k, e_k, 6)
        np.testing.assert_almost_equal(z, e_z, 6)
        np.testing.assert_almost_equal(p, e_p, 6)

if __name__ == '__main__':
    run_module_suite()
