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

# DEPRECATED MODULE

"""
Delsig based NTF synthesis
==========================

This modules provides code for the synthesis of the modulator NTF, based
on the routines in the DELSIG toolkit by R. Schreier.

.. deprecated:: 0.11.0
    Functions in ``delsig`` module are promoted to the ``NTFdesign`` module
    with slightly different names.
"""

from ..delsig import (synthesizeNTF, synthesizeChebyshevNTF,
                      clans)
