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
A collection of routines for the synthesis of Delta Sigma modulators.
=====================================================================

The delsig_optimize_NTF variable specifies whether to use optimization
techniques in the delsig_synthesizeNTF function.
"""

from . import filter_based
from . import weighting
from . import delsig
from . import psychoacoustic
from . import minmax

from .merit_factors import quantization_noise_gain
from .minmax import ntf_fir_minmax
from .delsig import (ntf_schreier, ntf_chebyshev, ntf_clans)
from .psychoacoustic import (ntf_dunn, ntf_fir_audio_weighting)

from . import tests as test_suite

from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
