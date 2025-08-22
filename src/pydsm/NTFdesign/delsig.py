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
# along with PyDSM.  If not, see <https://www.gnu.org/licenses/>.


"""
Delsig based NTF synthesis (:mod:`pydsm.NTFdesign.delsig`)
==========================================================

This modules provides code for the synthesis of the modulator NTF, based
on the routines in the DELSIG toolkit by R. Schreier.

These routines are also available in the :mod:`pydsm.delsig` and get
abbreviated names here.

.. currentmodule:: pydsm.NTFdesign.delsig


Functions
---------

.. autosummary::
   :toctree: generated/

    ntf_schreier    -- Design NTF with Schreier's method
    ntf_chebyshev   -- Design NTF based on Chebyshev type II form
    ntf_clans       -- Design NTF based on CLANS mehtod


Deprecated functions
--------------------

.. autosummary::
   :toctree: generated/

   synthesizeNTF           -- Alias for `ntf_schreier`
   synthesizeChebyshevNTF  -- Alias for `ntf_chebyshev`
   clans                   -- Alias for `ntf_clans`
"""

from __future__ import division, print_function

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
    """
    Alias of :func:`ntf_schreier`

    .. deprecated:: 0.11.0
        Function has been moved to the :mod:`pydsm.NTFdesign` module with
        name :func:`ntf_schreier`.
    """
    warn("Function superseded by ntf_schreier in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_schreier(order, osr, opt, H_inf, f0, **options)


def synthesizeChebyshevNTF(order=3, osr=64, opt=0, H_inf=1.5, f0=0.0):
    """
    Alias of :func:`ntf_chebyshev`

    .. deprecated:: 0.11.0
        Function has been moved to the :mod:`pydsm.NTFdesign` module with
        name :func:`ntf_chebyshev`.
    """
    warn("Function superseded by ntf_chebyshev in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_chebyshev(order, osr, opt, H_inf, f0)


def clans(order=4, osr=64, nq=5, rmax=0.95, opt=0, **options):
    """
    Alias of :func:`ntf_clans`

    .. deprecated:: 0.11.0
        Function has been moved to the :mod:`pydsm.NTFdesign` module with
        name :func:`ntf_clans`.
    """
    warn("Function superseded by ntf_clans in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_clans(order, osr, nq, rmax, opt, **options)
