# -*- coding: utf-8 -*-

# Copyright (c) 2014, Sergio Callegari
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

from ...correlations import raw_acorr

__all__ = ["q0_from_filter_ir"]


def q0_from_filter_ir(P, ir):
    """
    Legacy version of :func:`NTFDesign.q0_weighting`

    Legacy method based on correlation for the computation of the
    description matrix for the weighting function used in FIR NTF design.

    Parameters
    ----------
    P : int
        order of the FIR to be eventually synthesized
    ir : array_like
        output filter description impulse response

    Returns
    -------
    q0 : ndarray
        the first row of the matrix Q used in the NTF optimization

    Notes
    -----
    The Q matrix being synthesized has (P+1) times (P+1) entries.

    See Also
    --------
    scipy.integrate.quad : integrator used internally.
        For the meaning of the integrator parameters.
    """
    # Manage optional parameters
    return raw_acorr(ir, P)
