# -*- coding: utf-8 -*-

# Copyright (c) 2014, Sergio Callegari
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

from pydsm.NTFdesign.helpers import spread_fir_uc_zeros
from pydsm.NTFdesign.merit_factors import quantization_noise_gain
from pydsm.relab import cplxpair
from pydsm.delsig import ds_optzeros


class TestNTFdesignHelpers:
    def setUp(self):
        pass

    def test_uc_zeros(self):
        order = 7
        OSR = 32
        zeros1 = spread_fir_uc_zeros(order, OSR,
                                     quantization_noise_gain,
                                     cf_kwargs={'bounds': (0, 0.5/OSR)})
        zeros1 = cplxpair(zeros1)
        zeros2 = ds_optzeros(order)
        zeros2 = np.exp(1j*np.pi*zeros2/OSR)
        zeros2 = cplxpair(zeros2)
        np.testing.assert_almost_equal(zeros1, zeros2, 4)
