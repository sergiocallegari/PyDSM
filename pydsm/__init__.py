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
Library for the design and simulation of Delta-Sigma modulators
===============================================================

Copyright (c) 2012, Sergio Callegari
Some code ported from the DELSIG toolbox, copyrighted by R. Schreier.
"""

# Read version info
from ._version import __version__

# Promote some key functions to the pydsm namespace
from . import delsig
from . import NTFdesign
from . import simulation
from . import audio_weightings
from . import iso226

def test(level=1, verbosity=1):
    from numpy.testing import Tester
    return Tester().test(level, verbosity)