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
Acoustic weighting functions (:mod:`pydsm.audio_weightings`)
============================================================

Some standard acoustic weighting functions.

This module includes the A-, B- and C-weightings from the
ANSI Standards S1.4-1983 and S1.42-2001.

It also includes the D-weighting from the now withdrawn IEC 537.

It also includes the F-weighting proposed by R. A. Wannamaker.

The weighting functions can be expressed either in terms of
acoustic power or in terms of signal amplitude.

The weighting functions are also available in terms of a filter-based
implementation. In this case, be careful since no normalization is
present so that the gain at 1 kHz can be arbitrary. The filter
transfer function is referred to a signal amplitude weighting.

.. currentmodule:: pydsm.audio_weightings


Weighting functions
-------------------

.. autosummary::
   :toctree: generated/

    a_weighting
    b_weighting
    c_weighting
    d_weighting
    f_weighting


Filter implementation of weighting functions
--------------------------------------------

.. autodata:: a_zpk
   :annotation:
.. autodata:: b_zpk
   :annotation:
.. autodata:: c_zpk
   :annotation:
.. autodata:: d_zpk
   :annotation:
.. autodata:: f_zpk
   :annotation:


Normalization constants
-----------------------

.. autodata:: a_weighting_gain
   :annotation:
.. autodata:: b_weighting_gain
   :annotation:
.. autodata:: c_weighting_gain
   :annotation:
.. autodata:: d_weighting_gain
   :annotation:
.. autodata:: f_weighting_gain
   :annotation:


Notes
-----
The ANSI and IEC weightings are also described in Wikipedia [1]
and summarized in some illustrative web pages such as [2]_ and
[3]_. The F-weighting is documented in [4]_.

The filter-based implementation of the F-weighting is so high-order that
evaluation of the transfer function may require special care.

.. [1] Wikipedia (http://en.wikipedia.org/wiki/A-weighting)
.. [2] Cross spectrum (http://www.cross-spectrum.com/audio/weighting.html)
.. [3] Product Technology Parters "Noise Measurement Briefing"
   (http://www.ptpart.co.uk/noise-measurement-briefing/)
.. [4] Robert A. Wannamaker "Psychoacoustically Optimal Noise
   Shaping," J. Audio Eng. Soc., Vol. 40 No. 7/8 1992 July/August
"""

from __future__ import division, print_function

import numpy as np

__all__ = ["a_zpk", "a_weighting", "b_zpk", "b_weighting",
           "c_zpk", "c_weighting", "d_zpk", "d_weighting",
           "f_zpk", "f_weighting"]

a_zpk = (2*np.pi*np.asarray([0., 0., 0., 0.]),
         2*np.pi*np.asarray([-20.6, -20.6, -107.7, -739.9, -12200., -12200.]),
         (2*np.pi*12200.)**2)
"""A-weighting filter in zpk form."""

b_zpk = (2*np.pi*np.asarray([0., 0., 0.]),
         2*np.pi*np.asarray([-20.6, -20.6, -158.5, -12200., -12200.]),
         (2*np.pi*12200.)**2)
"""B-weighting filter in zpk form."""

c_zpk = (2*np.pi*np.asarray([0., 0.]),
         2*np.pi*np.asarray([-20.6, -20.6, -12200., -12200.]),
         (2*np.pi*12200.)**2)
"""C-weighting filter in zpk form."""

d_zpk = (2*np.pi*np.asarray([0., -519.8+876.2j, -519.8-876.2j]),
         2*np.pi*np.asarray([-282.7, -1160., -1712+2628j, -1712-2628j]),
         91104.32)
"""D-weighting filter in zpk form."""

f_zpk = (2*np.pi*np.asarray([0., 0., 0.,
                             -580+1030j, -580-1030j,
                             -3180+8750j, -3180-8750j,
                             -3180+8750j, -3180-8750j,
                             -3180+8750j, -3180-8750j]),
         2*np.pi*np.asarray([-180., -180., -180.,
                             -1630., -1630.,
                             -2510+3850j, -2510-3850j,
                             -2510+3850j, -2510-3850j,
                             -2510+3850j, -2510-3850j,
                             -2510+3850j, -2510-3850j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j,
                             -6620+14290j, -6620-14290j]),
         1.6810544531883432e+207)
"""F-weighting filter in zpk form."""

# Note: evaluating the transfer function of f_zpk may require special care
# since the high order implies that for many frequency values both the
# numerator and the denominator take very large values (in magnitude). Taking
# the ratio of large complex values may lead to overflow in numpy even if
# individually the numerator, the denominator and the result should not
# overflow.


def a_weighting(f, normal=True, power=True):
    """Returns the A-weighting as a function of frequency.

    Parameters
    ----------
    f : float or array of floats
        frequency where the weighting function is computed
    normal : bool
        whether the function should be normalized to have unit gain at
        1 kHz.
    power : bool
        whether the function should express the weighting in terms of
        acoustic power or signal amplitude

    Returns
    -------
    w : float or array of floats
        value of the weigting function
    """
    if power:
        return a_weighting(f, normal, power=False)**2
    w = (12200.0**2*f**4)/((f**2+20.6**2) *
                           np.sqrt((f**2+107.7**2) *
                           (f**2+737.9**2))*(f**2+12200.0**2))
    return w if not normal else w*a_weighting_gain

a_weighting_gain = 1/a_weighting(1000, normal=False, power=False)
"""Normalization gain to apply to A-weighting filter (namely, the
   attenuation of the filter at 1 kHz)"""


def b_weighting(f, normal=True, power=True):
    """Returns the B-weighting as a function of frequency.

    Parameters
    ----------
    f : float or array of floats
        frequency where the weighting function is computed
    normal : bool
        whether the function should be normalized to have unit gain at
        1 kHz.
    power : bool
        whether the function should express the weighting in terms of
        acoustic power or signal amplitude

    Returns
    -------
    w : float or array of floats
        value of the weigting function
    """
    if power:
        return b_weighting(f, normal, power=False)**2
    w = (12200.0**2*f**3)/((f**2+20.6**2) *
                           np.sqrt(f**2+158.5**2)*(f**2+12200.0**2))
    return w if not normal else w*b_weighting_gain

b_weighting_gain = 1/b_weighting(1000, normal=False, power=False)
"""Normalization gain to apply to B-weighting filter (namely, the
   attenuation of the filter at 1 kHz)"""


def c_weighting(f, normal=True, power=True):
    """Returns the C-weighting as a function of frequency.

    Parameters
    ----------
    f : float or array of floats
        frequency where the weighting function is computed
    normal : bool
        whether the function should be normalized to have unit gain at
        1 kHz.
    power : bool
        whether the function should express the weighting in terms of
        acoustic power or signal amplitude

    Returns
    -------
    w : float or array of floats
        value of the weigting function
    """
    if power:
        return c_weighting(f, normal, power=False)**2
    w = (12200.0**2*f**2)/((f**2+20.6**2)*(f**2+12200.0**2))
    return w if not normal else w*c_weighting_gain

c_weighting_gain = 1/c_weighting(1000, normal=False, power=False)
"""Normalization gain to apply to C-weighting filter (namely, the
   attenuation of the filter at 1 kHz)"""


def d_weighting(f, normal=True, power=True):
    """Returns the D-weighting as a function of frequency.

    Parameters
    ----------
    f : float or array of floats
        frequency where the weighting function is computed
    normal : bool
        whether the function should be normalized to have unit gain at
        1 kHz. This parameter is ignored, since this weighting function
        is always normalized.
    power : bool
        whether the function should express the weighting in terms of
        acoustic power or signal amplitude

    Returns
    -------
    w : float or array of floats
        value of the weigting function
    """
    if power:
        return d_weighting(f, normal, power=False)**2

    def h(f):
        return (((1037918.48-f**2)**2+1080768.16*f**2) /
                ((9837328.0-f**2)**2+11723776.0*f**2))

    return (f/6.8966888496476E-5 *
            np.sqrt(h(f)/((f**2+79919.29)*(f**2+1345600.0))))

d_weighting_gain = 1.
"""Normalization gain to apply to D-weighting filter (namely, the
   attenuation of the filter at 1 kHz)"""


def f_weighting(f, normal=True, power=True):
    """Returns the F-weighting as a function of frequency.

    Parameters
    ----------
    f : float or array of floats
        frequency where the weighting function is computed
    normal : bool
        whether the function should be normalized to have unit gain at
        1 kHz.
    power : bool
        whether the function should express the weighting in terms of
        acoustic power or signal amplitude

    Returns
    -------
    w : float or array of floats
        value of the weigting function

    Notes
    -----
    The F-weighting function is documented in [1]_.

    References
    ----------
    .. [1] Robert A. Wannamaker "Psychoacoustically Optimal Noise Shaping,"
       J. Audio Eng. Soc., Vol. 40 No. 7/8 1992 July/August
    """
    if not power:
        return np.sqrt(f_weighting(f, normal, power=True))
    fx = f/1000.
    g = 2.536e-5
    z1 = fx**2
    z2 = ((0.58**2)+(1.03**2)-z1)**2 + 4.0*(0.58**2)*z1
    z3 = ((3.18**2)+(8.75**2)-z1)**2 + 4.0*(3.18**2)*z1
    p1 = 0.18**2+z1
    p2 = 1.63**2+z1
    p3 = ((2.51**2)+(3.85**2)-z1)**2 + 4.0*(2.51**2)*z1
    p4 = ((6.62**2)+(14.29**2)-z1)**2 + 4.0*(6.62**2)*z1
    w = ((g*((z1**3)*z2*(z3**3)) /
         ((p1**3)*(p2**2)*(p3**4))*((1e5/p4)**20)))
    return w if not normal else w*f_weighting_gain

# Set normalization gain
f_weighting_gain = 1/f_weighting(1000, normal=False, power=True)
"""Normalization gain to apply to F-weighting filter (namely, the
   attenuation of the filter at 1 kHz)"""
