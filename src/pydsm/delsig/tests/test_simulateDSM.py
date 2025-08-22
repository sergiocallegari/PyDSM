# -*- coding: utf-8 -*-

# Copyright (c) 2012–2024, Sergio Callegari
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
# along with PyDSM.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import division, print_function

import numpy as np
import importlib_resources
from pydsm.delsig import simulateDSM

__all__ = ["TestSimulateDSM"]


class TestSimulateDSM:

    def setUp(self):
        pass

    def test_default(self):
        with (importlib_resources.files('pydsm.delsig')
              .joinpath('tests/Data/test_simulateDSM_0.npz')
              .open('rb')) as f:
            d = np.load(f)['arr_0']
        # Take H as in H = synthesizeNTF(5, 32, 1)
        H = (np.array([0.99604531+0.08884669j,  0.99604531-0.08884669j,
                       0.99860302+0.05283948j,  0.99860302-0.05283948j,
                       1.00000000+0.j]),
             np.array([0.80655696+0.11982271j,  0.80655696-0.11982271j,
                       0.89807098+0.21981939j,  0.89807098-0.21981939j,
                       0.77776708+0.j]),
             1)
        N = 8192
        f = 85
        u = 0.5*np.sin(2.*np.pi*f/N*np.arange(N))
        v, d1, d2, d3 = simulateDSM(u, H)
        np.testing.assert_equal(v, d)
