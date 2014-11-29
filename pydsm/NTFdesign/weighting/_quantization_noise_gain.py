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


# Following part is deprecated

from __future__ import division, print_function

from ..merit_factors import quantization_noise_gain
from warnings import warn
from ...exceptions import PyDsmDeprecationWarning

__all__ = ["quantization_weighted_noise_gain"]


def quantization_weighted_noise_gain(NTF, w=None, bounds=(0, 0.5),
                                     **options):
    warn("Function superseded by quantization_noise_gain in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return quantization_noise_gain(NTF, w, bounds, **options)

quantization_weighted_noise_gain.__doc__ = \
    quantization_noise_gain.__doc__ + """
    .. deprecated:: 0.11.0
        Function has been moved to the ``NTFdesign`` module with name
        ``quantization_noise_gain``.
    """

quantization_weighted_noise_gain.default_options = \
    quantization_noise_gain.default_options
