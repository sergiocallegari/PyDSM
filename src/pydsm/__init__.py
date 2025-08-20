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

"""PyDSM Toolbox for the design and simulation of ΔΣ modulators
============================================================

PyDSM © 2012–2025 Sergio Callegari. All rights reserved.

Some code portions confined in the :mod:`pydsm.delsig` module are
ported or reimplemented from the DELSIG_ toolbox for MATLAB,
copyrighted by R. Schreier.

PyDSM focuses on the design of Noise Transfer Functions for Digital ΔΣ
Modulators and in the simulation of the resulting modulators.


Main modules
------------

.. toctree::
   :maxdepth: 1

   pydsm.delsig
   pydsm.NTFdesign
   pydsm.simulation
   pydsm.audio_weightings
   pydsm.iso226

Utility modules
---------------

.. toctree::
   :maxdepth: 1

   pydsm.correlations
   pydsm.ft
   pydsm.ir
   pydsm.relab
   pydsm.utilities


Error handling
--------------

.. toctree::
   :maxdepth: 1

   pydsm.exceptions

.. _DELSIG:
   http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox

"""

# Read version info
from ._version import __version__

# Promote some key functions to the pydsm namespace
from . import delsig
from . import NTFdesign
from . import simulation
from . import audio_weightings
from . import iso226

from ._pytesttester import PytestTester
test = PytestTester(__name__)
bench = PytestTester(__name__, mode='benchmarks')
del PytestTester
