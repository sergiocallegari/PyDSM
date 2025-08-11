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

"""
Simulation of Delta Sigma modulators (:mod:`pydsm.simulation`)
==============================================================

This module provides some functions for the simulation of of ΔΣ
modulators. Currently, the functions included in this module are just
alternative entry points for functions in the :mod:`pydsm.delsig` module.

.. currentmodule:: pydsm.simulation


Functions
---------

.. autosummary::
   :toctree: generated/

   simulateDSM   -- Delta sigma modulator simulation
   ds_quantize   -- quantization function
"""

# Promote some functions/global variables to the simulation namespace
from .delsig import simulateDSM
from .delsig import ds_quantize

__all__ = ['simulateDSM', 'ds_quantize']
