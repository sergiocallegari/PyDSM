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

u"""
Output filter based NTF synthesis
=================================

This modules provides code for the synthesis of the modulator NTF, based
on the filter placed after the modulator for the removal of the quantization
noise.

The algorithms used in these routines are extensively described in

Sergio Callegari, Federico Bizzarri "Output Filter Aware Optimization of the
Noise Shaping Properties of ΔΣ Modulators via Semi-Definite Programming",
IEEE Transactions on Circuits and Systems I: Regular Papers.
"""

from ._q0_from_filter import *
from ._quantization_noise_gain import *
from ._ntf_from_filter import *