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

u"""
Synthesis of the NTF of Delta Sigma modulators (:mod:`pydsm.NTFdesign`)
=======================================================================

This module provides some strategies for the design of the Noise
Transfer Function of ΔΣ modulators. There are both functions that are
specific to PyDSM and entry points to functions in the delsig module
of PyDSM (:mod:`pydsm.delsig`).

.. currentmodule:: pydsm.NTFdesign


Key functions
-------------

.. function:: ntf_schreier()

   shorthand for :func:`delsig.ntf_schreier`

.. function:: ntf_chebyshev()

   shorthand for :func:`delsig.ntf_chebyshev`

.. function:: ntf_clans()

   shorthand for :func:`delsig.ntf_clans`

.. function:: ntf_fir_weighting()

   shorthand for :func:`weighting.ntf_fir_weighting`

.. function:: ntf_hybrid_weighting()

   shorthand for :func:`weighting.ntf_hybrid_weighting`

.. function:: ntf_fir_minmax()

   shorthand for :func:`minmax.ntf_fir_minmax`

.. function:: ntf_dunn()

   shorthand for :func:`psychoacoustic.ntf_dunn`

.. function:: ntf_fir_audio_weighting()

   shorthand for :func:`psychoacoustic.ntf_fir_audio_weighting`

.. function:: mult_weightings()

   shorthand for :func:`weighting.mult_weightings`

.. function:: quantization_noise_gain()

   shorthand for :func:`merit_factors.quantization_noise_gain`


Submodules
----------

:mod:`pydsm.NTFdesign.delsig`
  NTF synthesis functions equivalent to those in :mod:`pydsm.delsig`.

:mod:`pydsm.NTFdesign.weighting`
  NTF synthesis techniques that take as their input either a weighting
  function (indicating the cost of quantization noise power versus frequency)
  or a specification of the filter in charge of removing the quantization
  noise.

:mod:`pydsm.NTFdesign.minmax`
  NTF synthesis techniques based on a minmax approach.

:mod:`pydsm.NTFdesign.psychoacoustic`
  NTF synthesis techniques for audio modulators that result in a noise shaping
  that take into account psychoacoustics.

:mod:`pydsm.NTFdesign.merit_factors`
  Functions for determining merit factors about NTFs.

:mod:`pydsm.NTFdesign.helpers`
  Helper functions


Legacy submodule
----------------

:mod:`pydsm.NTFdesign.legacy`
  Functions that are now superseded but that are worth keeping around for
  reproducing published results


Deprecated submodules
---------------------

:mod:`pydsm.NTFdesign.filter_based`
  Alternate entry points for some functions.
"""


from .merit_factors import quantization_noise_gain
from .minmax import ntf_fir_minmax
from .delsig import ntf_schreier, ntf_chebyshev, ntf_clans
from .psychoacoustic import ntf_dunn, ntf_fir_audio_weighting
from .weighting import (ntf_fir_weighting, ntf_hybrid_weighting,
                        mult_weightings)

from .._pytesttester import PytestTester
test = PytestTester(__name__)
bench = PytestTester(__name__, mode='benchmarks')
del PytestTester
