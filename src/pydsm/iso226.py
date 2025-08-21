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
# along with PyDSM.  If not, see <https://www.gnu.org/licenses/>.

"""
ISO 226 loudness countours (:mod:`pydsm.iso226`)
================================================

Loudness contours from ISO 226.

Contours are returned both as tables of data and as contour functions

.. currentmodule:: pydsm.iso226


Functions returning ISO 226 contours
------------------------------------

.. autosummary::
   :toctree: generated/

   iso226_spl_contour -- Equal loudness contour (tabled)
   iso226_spl_itpl  -- Interpolated equal loudness contour


Functions computing loudness/acoustic pressure
----------------------------------------------

.. autosummary::
   :toctree: generated/

   tabled_L_p  -- Return table of sound pressure levels for given loudness
   tabled_L_N  -- Return table of perceived loudness for given sound pressure


Functions returning data tabled in the standard
-----------------------------------------------

.. autosummary::
   :toctree: generated/

   tabled_f  -- Return table of frequencies in ISO 226
   tabled_alpha_f  -- Return table of exponents for loudness perception
   tabled_L_U  -- Return table of magnitudes of the linear transfer function
   tabled_T_f  -- Return table of thresholds of hearing


Notes
-----
This module uses data from the latest revision of ISO 226 [1]_.
For reference, also consider [2]_.

The ISO standard provides the equal loudness contours as tabled data.
Tables end at 12.5 kHz. Above this frequency equal-loudness-level data are
relatively scarce and tend to be variable [3]_. Yet, it is known
that the human ear has a precipitous decline in sensitivity with increasing
frequency above 15 kHz, to the point that at about 20 kHz the percieved sound
becomes negligible (> 100 dB attenuation) [4]_. For this reason,
this module includes the possibility of delivering some modified ISO contours
where the tabled data is augmented by creating a new data point at 20 kHz
where the behavior found at 20 Hz is replicated.

.. [1] ISO 226:2003 "Acoustics - Normal equal-loudness-level contours"
.. [2] Jeff Tackett, "ISO 226 Equal-Loudness-Level Contour Signal,"
   2005
   (https://www.mathworks.com/matlabcentral/fileexchange/7028)
.. [3] Yoiti Suzuki et al, "Precise and Full-range Determination of
   Two-dimensional Equal Loudness Contours," 2003
   (http://www.mp3-tech.org/programmer/docs/IS-01Y-E.pdf)
.. [4] Robert A. Wannamaker "Psychoacoustically Optimal Noise
   Shaping", J. Audio Eng. Soc., Vol. 40, N. 7/8, 1992 July/August
"""

from __future__ import division, print_function

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

__all__ = ["tabled_f", "tabled_alpha_f", "tabled_L_U", "tabled_T_f",
           "tabled_L_p", "tabled_L_N",
           "iso226_spl_contour", "iso226_spl_itpl"]

# Tabled ISO 226 parameters
tbl_f = np.asarray(
    [20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400,
     500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300,
     8000, 10000, 12500])
tbl_alpha_f = np.asarray(
    [0.532, 0.506, 0.480, 0.455, 0.432, 0.409, 0.387, 0.367, 0.349, 0.330,
     0.315, 0.301, 0.288, 0.276, 0.267, 0.259, 0.253, 0.250, 0.246, 0.244,
     0.243, 0.243, 0.243, 0.242, 0.242, 0.245, 0.254, 0.271, 0.301])
tbl_L_U = np.asarray(
    [-31.6, -27.2, -23.0, -19.1, -15.9, -13.0, -10.3, -8.1, -6.2, -4.5,
     -3.1, -2.0, -1.1, -0.4, 0.0, 0.3, 0.5, 0.0, -2.7, -4.1, -1.0, 1.7,
     2.5, 1.2, -2.1, -7.1, -11.2, -10.7, -3.1])
tbl_T_f = np.asarray(
    [78.5, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1, 17.9, 14.4,
     11.4, 8.6, 6.2, 4.4, 3.0, 2.2, 2.4, 3.5, 1.7, -1.3, -4.2, -6.0, -5.4,
     -1.5, 6.0, 12.6, 13.9, 12.3])


def tabled_f(hfe=False):
    """Table of frequencies in ISO 226.

    Parameters
    ----------
    hfe : bool
        whether the table should be augmented with a data point
        at 20 kHz (High-Frequency-Enhanced table)

    Returns
    -------
    f : array of floats
        the frequency table.
    """
    return np.append(tbl_f, 20E3) if hfe else tbl_f


def tabled_alpha_f(hfe=False):
    """Table of exponents for loudness perception in ISO 226.

    Parameters
    ----------
    hfe : bool
        whether the table should be augmented with a data point
        at 20 kHz (High-Frequency-Enhanced table)

    Returns
    -------
    alpha_f : array of floats
        the exponents table.
    """
    return np.append(tbl_alpha_f, tbl_alpha_f[0]) if hfe else tbl_alpha_f


def tabled_L_U(hfe=False):
    """Table of magnitudes of the linear transfer function in ISO 226.

    Parameters
    ----------
    hfe : bool
        whether the table should be augmented with a data point
        at 20 kHz (High-Frequency-Enhanced table)

    Returns
    -------
    L_U : array of floats
        the magnitudes table.

    Notes
    -----
    The returned values are the magnitude of the linear transfer function
    normalized at 1 kHz.
    """
    return np.append(tbl_L_U, tbl_L_U[0]) if hfe else tbl_L_U


def tabled_T_f(hfe=False):
    """Table of thresholds of hearing in ISO 226.

    Parameters
    ----------
    hfe : bool
        whether the table should be augmented with a data point
        at 20 kHz (High-Frequency-Enhanced table)

    Returns
    -------
    T_f : array of floats
        the thresholds table.
    """
    return np.append(tbl_T_f, tbl_T_f[0]) if hfe else tbl_T_f


# Check that it works fine when L_N is array
def tabled_A_f(L_N, hfe=False):
    """Table of A_f values for given loundess in ISO 226.

    Parameters
    ----------
    L_N : float
        percieved loudness level in phons
    hfe : bool
        whether the table should be augmented with a data point
        at 20 kHz (High-Frequency-Enhanced table)

    Returns
    -------
    A_f : array of floats
        the A_f table.

    Notes
    -----
    1 phon is 1 dB_SPL (sound pressure level) at 1 kHz. Sound pressure levels
    are measured in dBs by referring to a reference pressure level P0 (close
    to the hearing threshold at 1 kHz and set to 20 uPa RMS).
    """
    A_f = (4.47E-3*(10.0**(0.025*L_N)-1.15) +
           (0.4*10.0**((tbl_T_f+tbl_L_U)/10.0-9.0))**tbl_alpha_f)
    return np.append(A_f, A_f[0]) if hfe else A_f


# Check that it works fine when L_N is array
def tabled_L_p(L_N, hfe=False):
    """Table of sound pressure levels for given loudness in ISO 226.

    This function returns a table according to ISO 226 sect 4.1.

    Parameters
    ----------
    L_N : float
        percieved loudness level in phons
    hfe : bool
        whether the table should be augmented with a data point
        at 20 kHz (High-Frequency-Enhanced table)

    Returns
    -------
    L_p : array of floats
        the sound pressure level table. Sound pressure levels are returned
        in DB_SPL

    Notes
    -----
    1 phon is 1 dB_SPL (sound pressure level) at 1 kHz. Sound pressure levels
    are measured in dBs by referring to a reference pressure level P0 (close
    to the hearing threshold at 1 kHz and set to 20 uPa RMS).
    """
    L_p = (10.0/tbl_alpha_f)*np.log10(tabled_A_f(L_N))-tbl_L_U + 94.0
    return np.append(L_p, L_p[0]) if hfe else L_p


# Check that it works fine when L_P is array
def tabled_B_f(L_p, hfe=False):
    """Table of B_f values for given sound pressure in ISO 226.

    Parameters
    ----------
    L_p : float
        sound pressure level in dB_SPL
    hfe : bool
        whether the table should be augmented with a data point
        at 20 kHz (High-Frequency-Enhanced table)

    Returns
    -------
    B_f : array of floats
        the B_f table.

    Notes
    -----
    Sound pressure levels are measured in dBs by referring to a
    reference pressure level P0 (close to the hearing threshold at 1 kHz and
    set to 20 uPa RMS).
    """
    B_f = ((0.4*10**(L_p+tbl_L_U)/10.-9.)**tbl_alpha_f -
           (0.4*10**(tbl_T_f+tbl_L_U)/10.-9.)**tbl_alpha_f +
           0.005135)
    return np.append(B_f, B_f[0]) if hfe else B_f


# Check that it works fine when L_N is array
def tabled_L_N(L_p, hfe=False):
    """Table of perceived loudness levels for given sound pressure in ISO 226.

    This function returns a table according to ISO 226 sect 4.2.

    Parameters
    ----------
    L_p : float
        sound pressure level in dB_SPL
    hfe : bool
        whether the table should be augmented with a data point
        at 20 kHz (High-Frequency-Enhanced table)

    Returns
    -------
    L_N : array of floats
        the perceived loudness level table. Loudness levels are returned
        in phons

    Notes
    -----
    1 phon is 1 dB_SPL (sound pressure level) at 1 kHz. Sound pressure levels
    are measured in dBs by referring to a reference pressure level P0 (close
    to the hearing threshold at 1 kHz and set to 20 uPa RMS).
    """
    L_N = 40*np.log10(tabled_B_f(L_p))+94.0
    return np.append(L_N, L_N[0]) if hfe else L_N


def iso226_spl_contour(L_N=40, hfe=False):
    """Generates an equal loudness contour as described in ISO 226.

    This function returns the control points describing the equal loudness
    contour for the input phon level, according to ISO 226 sect 4.1.

    Parameters
    ----------
    L_N : float, optional
        perceived loudness level in phons.
    hfe : bool
        whether the table should be augmented with a data point
        at 20 kHz (High-Frequency-Enhanced table)

    Returns
    -------
    f : ndarray
        frequencies where the loudness is evaluated.
    spl : ndarray
        equivalent sound pressure level at the frequencies f.

    Notes
    -----
    1 phon is 1 dB_SPL (sound pressure level) at 1 kHz. Sound pressure levels
    are measured in dBs by referring to a reference pressure level P0 (close
    to the hearing threshold at 1 kHz and set to 20 uPa RMS).

    The valid input phon range is 0-90 dB_SPL. Above 80 dB, only the
    frequency range 20-4000 Hz is significant.
    """
    # Check for valid input range
    if L_N < 0 or L_N > 90:
        raise ValueError('Parameter L_N out of bounds [0-90].')
    # Derive sound pressure level from loudness level (ISO 226 sect 4.1)
    return tabled_f(hfe), tabled_L_p(L_N, hfe)


def iso226_spl_itpl(L_N=40, hfe=False, k=3):
    """Generates an interpolation of an equal loudness contour.

    This function returns an interpolation object describing the equal
    loudness contour for the input phon level, according to ISO 226 sect 4.1.

    Parameters
    ----------
    L_N : float, optional
        perceived loudness level in phons.
    hfe : bool
        whether the table should be augmented with a data point
        at 20 kHz (High-Frequency-Enhanced table)
    k : int
        interpolation order

    Returns
    -------
    itpl : univariate interpolation object
        function-like object that takes a frequency f as its input and returns
        the equivalent sound pressure level at f

    Notes
    -----
    1 phon is 1 dB_SPL (sound pressure level) at 1 kHz. Sound pressure levels
    are measured in dBs by referring to a reference pressure level P0 (close
    to the hearing threshold at 1 kHz and set to 20 uPa RMS).

    The valid input phon range is 0-90 dB_SPL. Above 80 dB, only the
    frequency range 20-4000 Hz is significant.
    """
    ff, yy = iso226_spl_contour(L_N, hfe)
    return InterpolatedUnivariateSpline(ff, yy, k=k)
