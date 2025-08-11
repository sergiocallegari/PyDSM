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
Noise weighting NTF synthesis (:mod:`pydsm.NTFdesign.weighting`)
================================================================

This modules provides code for the synthesis of the modulator NTF, based
on a noise weighting function.

.. currentmodule:: pydsm.NTFdesign.weighting


Functions
---------

.. autosummary::
   :toctree: generated/

    mult_weightings      -- Product of weighting functions
    q0_weighting         -- Q matrix from a noise weighting function
    ntf_fir_from_q0      -- design FIR NTF from matrix Q
    ntf_hybrid_weighting -- design hybrid NTF from weighting function
    ntf_fir_weighting    -- design FIR NTF from weighting function


Deprecated Functions
--------------------

.. autosummary::
   :toctree: generated/

   quantization_weighted_noise_gain    -- Alias of `quantization_noise_gain`
   q0_from_noise_weighting             -- Alias of `q0_weighting`
   synthesize_ntf_from_q0              -- Alias of `ntf_fir_from_q0`
   synthesize_ntf_from_noise_weighting -- Alias of `ntf_fir_weighting`
"""

from ._fir_weighting import *
from ._quantization_noise_gain import *

__all__ = (_fir_weighting.__all__ + _quantization_noise_gain.__all__
           + ['test'])

from ..._pytesttester import PytestTester
test = PytestTester(__name__)
del PytestTester
