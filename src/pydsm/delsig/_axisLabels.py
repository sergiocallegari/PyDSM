# -*- coding: utf-8 -*-

# Copyright (c) 2015, Sergio Callegari
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

# This file includes code ported from the DELSIG Matlab toolbox
# (see http://www.mathworks.com/matlabcentral/fileexchange/19)
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

from __future__ import division

import numpy as np

__all__ = ["axisLabels"]


def axisLabels(points, incr):
    """
    Generate alphanumeric axis labels.

    If the parameter `incr` is an integer, the function returns an array
    of strings corresponding to the values specified in `points`,
    stepping through them with indexes 0, incr, 2*incr, ...

    If the parameter `incr` is a couple of integers `(start, incr)`,
    the stepping occurs with indexes start, start+incr, start+2*incr, ...

    Parameters
    ----------
    points : array like of floats
        the axis points to label
    incr : int or couple of ints
        the increment to use in labelling

    Returns
    -------
    s : list of strings
        list of labels

    Notes
    -----

    In the original DELSIG function, the parameter `points` is named `range`.

    Values in `points` less than 1e-6 are truncated to zero. This value is
    hardwired in the function.
    """
    points = np.asarray(points).flatten()
    # Flatten assures that a copy is taken, so that the following assignement
    # does not alter the original vector
    points[np.abs(points) < 1e-6] = 0
    if isinstance(incr, int):
        start = 0
        step = incr
    elif (hasattr(incr, '__len__') and
          len(incr) == 2 and
          np.all([isinstance(x, int) for x in incr])):
        start = incr[0]
        step = incr[1]
    else:
        raise TypeError('Parameter incr must be integer or couple of integer '
                        'values')
    return ['%g' % points[i] for i in range(start, len(points), step)]
