# -*- coding: utf-8 -*-

# Copyright (c) 2012–2025, Sergio Callegari
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

# This file includes code ported from the DELSIG Matlab toolbox
# (see https://www.mathworks.com/matlabcentral/fileexchange/19)
# covered by the following copyright and permission notice
#
# Copyright (c) 2009 Richard Schreier
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the distribution
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""
Code ported/reimplemented from the DELSIG toolbox (:mod:`pydsm.delsig`)
=======================================================================

The DELSIG toolbox for MATLAB is copyrighted by R. Schreier and
distributed under the BSD license.

.. currentmodule:: pydsm.delsig


Key functions
-------------

.. autosummary::
   :toctree: generated/

   synthesizeNTF
   clans
   synthesizeChebyshevNTF
   simulateDSM

Other selected functions
------------------------

Delta sigma utilities
.....................

.. autosummary::
   :toctree: generated/

   partitionABCD
   rmsGain

General utilities
.................

.. autosummary::
   :toctree: generated/

   dbv
   dbp
   undbv
   undbp
   dbm
   undbm
   rms

Graphing
........

.. autosummary::
   :toctree: generated/

   plotPZ
   axisLabels


Plumbing
--------

.. autosummary::
   :toctree: generated/

   ds_synNTFobj1
   ds_f1f2
   ds_optzeros
   dsclansNTF
   padl
   padr
   padt
   padb
   evalTF
   evalRPoly

"""

__delsig_version__ = "7.4"

# The delsig module reflects the flat organization of the original DELSIG
from ._axisLabels import *
from ._decibel import *
from ._ds import *
from ._padding import *
from ._tf import *
from ._plot import *
from ._synthesizeNTF import *
from ._synthesizeChebyshevNTF import *
from ._clans import *
from ._dsclansNTF import *
from ._simulateDSM import *
from ._simulateDSM_scipy import *
from ._partitionABCD import *
from ._rmsGain import *
from ._rms import *

from .._pytesttester import PytestTester
test = PytestTester(__name__)
bench = PytestTester(__name__, mode='benchmarks')
del PytestTester
