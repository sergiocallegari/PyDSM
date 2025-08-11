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
# along with PyDSM.  If not, see <http://www.gnu.org/licenses/>.

u"""
Legacy NTF synthesis methods (:mod:`pydsm.NTFDesign.legacy`)
============================================================

This modules provides code for the synthesis of the modulator NTF, based
on methods that are now superseded but that are worth keeping around for
reproducing published results.

.. currentmodule:: pydsm.NTFdesign.legacy


Functions
---------

.. autosummary::
   :toctree: generated/

   q0_from_filter_ir                -- Legacy `q0_weighting`
   quantization_noise_gain_by_conv  -- Legacy `quantization_noise_gain`
"""

from ._quantization_noise_gain import *
from ._q0_from_ir import *

__all__ = (_quantization_noise_gain.__all__ + _q0_from_ir.__all__
           + ['test'])

from ..._pytesttester import PytestTester
test = PytestTester(__name__)
del PytestTester
