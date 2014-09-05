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

"""
Functions to obtain the matrix used in the NTF optimization
===========================================================

These routines return the matrix from a noise weigthing function.

The matrix is symmetric Toeplitz, thus actually described by its first
row only, which is what the routines here return.
"""

from ...ft import idtft_hermitian
import numpy as np

__all__=["q0_from_noise_weighting"]

def q0_from_noise_weighting(P, noise_weighting, \
    integrator_params={'epsabs':1E-14, 'epsrel':1E-9}):
    """
    Computes Q matrix from a noise weighting function

    Parameters
    ----------
    P : int
        order of the FIR to be eventually synthesized
    noise_weighting : callable
        function of f acting as a noise weighting.
        f is normalized between 0 and 0.5

    Returns
    -------
    q0 : ndarray
        the first row of the matrix Q used in the NTF optimization

    Other parameters
    ----------------
    integrator_params : dict, optional
        the controlling parameters for the numerical integrator
        (see `scipy.integrate.quad`)

    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.
    """
    ac=lambda t:idtft_hermitian(noise_weighting, t, integrator_params)
    return np.asarray(map(ac,np.arange(P+1)))
