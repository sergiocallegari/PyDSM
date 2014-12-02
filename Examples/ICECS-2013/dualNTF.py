# -*- coding: utf-8 -*-

# Copyright (c) 2013, Sergio Callegari
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

from __future__ import division, print_function

import numpy as np
from scipy.signal import tf2zpk, zpk2tf

__all__ = ['mirroredNTF']


def mirroredNTF(ntf):
    """Creates a symmetric noise transfer function from a prototype.

    This function takes a prototype single band NTF for a LP modulator
    and returns a symmetric NTF with signal bands at the bottom and at
    the top of the available frequency ranges.

    Note that the way in which the mirroring happens make the resulting
    NTF peak at about the squared value as the peak of the original
    prototype. Hence the prototype should be designed with an H_inf value
    that is the square root of the desired one.

    Parameters
    ----------
    ntf : tuple
        input noise tranfer function in zpk form

    Returns
    -------
    ntf2 : tuple
        output noise transfer function in zpk form
    """
    zpk_mode = True
    if len(ntf) != 3:
        zpk_mode = False
        ntf = tf2zpk(*ntf)
    # Take the opposite of poles and zeros to get an HP filter
    ntf_flipped = (-ntf[0], -ntf[1], ntf[2])
    # Take the product of the two
    ntf_mirrored = (np.hstack((ntf[0], ntf_flipped[0])),
                    np.hstack((ntf[1], ntf_flipped[1])), 1)
    if not zpk_mode:
        ntf_mirrored = zpk2tf(*ntf_mirrored)
    return ntf_mirrored
