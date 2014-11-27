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

# This file includes code ported from the "Frequency-Domain Min-Max
# Optimization for Delta-Sigma Modulators" Matlab toolbox by Masaaki
# Nagahara (see http://www.mathworks.com/matlabcentral/fileexchange/36187)
# covered by the following copyright and permission notice
#
# Copyright (c) 2012, Masaaki Nagahara
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the distribution
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# DEPRECATED MODULE

from warnings import warn
from ..exceptions import PyDsmDeprecationWarning
from . import _minmax

__doc__ = _minmax.__doc__ + """
    .. deprecated:: 0.11.0
    Functions in ``minmax`` module are promoted to the ``NTFdesign`` module
    with slightly different names.
    """

__all__ = ['synthesize_ntf_minmax']


def synthesize_ntf_minmax(order=32, osr=32, H_inf=1.5, f0=0, zf=False,
                          **options):
    warn("Function superseded by ntf_fir_minmax in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return _minmax.ntf_fir_minmax(order, osr, H_inf, f0, zf, **options)

synthesize_ntf_minmax.__doc__ = _minmax.ntf_fir_minmax.__doc__ + """
    .. deprecated:: 0.11.0
    Function has been moved to the ``NTFdesign`` module with name
    ``ntf_fir_minmax``.
    """

synthesize_ntf_minmax.default_options = \
    _minmax.ntf_fir_minmax.default_options.copy()
