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
Delsig based NTF synthesis
==========================

This modules provides code for the synthesis of the modulator NTF, based
on the routines in the DELSIG toolkit by R. Schreier.

The routines have abbreviated names here.
"""

from __future__ import division, print_function

import numpy as np
from ..delsig import (synthesizeNTF as ntf_schreier,
                      synthesizeChebyshevNTF as ntf_chebyshev,
                      clans as ntf_clans)
from warnings import warn
from ..exceptions import PyDsmDeprecationWarning

__all__ = ['ntf_schreier', 'ntf_chebyshev', 'ntf_clans',
           'synthesizeNTF', 'synthesizeChebyshevNTF', 'clans']


# Following part is deprecated

def synthesizeNTF(order=3, osr=64, opt=0, H_inf=1.5, f0=0.0,
                  **options):
    warn("Function superseded by ntf_schreier in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_schreier(order, osr, opt, H_inf, f0, **options)

synthesizeNTF.default_options = ntf_schreier.default_options

synthesizeNTF.__doc__ = ntf_schreier.__doc__ + """
    .. deprecated:: 0.11.0
        Function has been moved to the ``NTFdesign`` module with
        name ``ntf_schreier``.
    """


def synthesizeChebyshevNTF(order=3, osr=64, opt=0, H_inf=1.5, f0=0.0):
    warn("Function superseded by ntf_chebyshev in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_chebyshev(order, osr, opt, H_inf, f0)

synthesizeChebyshevNTF.__doc__ = ntf_chebyshev.__doc__ + """
    .. deprecated:: 0.11.0
        Function has been moved to the ``NTFdesign`` module with
        name ``ntf_chebyshev``.
    """


def clans(order=4, osr=64, nq=5, rmax=0.95, opt=0, **options):
    warn("Function superseded by ntf_clans in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_clans(order, osr, nq, rmax, opt, **options)

clans.default_options = ntf_clans.default_options

clans.__doc__ = ntf_clans.__doc__ + """
    .. deprecated:: 0.11.0
        Function has been moved to the ``NTFdesign`` module with
        name ``ntf_clans``.
    """
