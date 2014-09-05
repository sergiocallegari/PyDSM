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
Computation of the noise power gain through the NTF and the output filter
=========================================================================
"""

import numpy as np
from scipy.integrate import quad
from ...delsig import evalTF

__all__ = ["quantization_weighted_noise_gain"]

def quantization_weighted_noise_gain(NTF, w):
    r"""
    Computes the quantization noise power gain

    Parameters
    ----------
    NTF : tuple
        NTF definition in zpk or nd form
    w : callable with argument f in [0,1/2]
        noise weighting function

    Returns
    -------
    a : real
        noise power gain

    Notes
    -----
    The computation is practiced as

    .. math::
        2\int_{f=0}^{\frac{1}{2}}
        \left|\mathit{NTF}
        \left(\mathrm{e}^{\mathrm{i} 2\pi f}\right)\right|^2 w(f) df
    """
    return 2*quad(lambda f: \
        np.abs(evalTF(NTF,np.exp(2j*np.pi*f)))**2 * w(f), 0, 0.5)[0]
