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
Fourier transform related routines
==================================

Functions to compute the fft and the dtft.
"""

from __future__ import division, print_function

import numpy as np
import scipy as sp
__import__("scipy.fftpack")
__import__("scipy.integrate")

__all__ = ["fft_centered", "dtft", "dtft_hermitian", "idtft",
           "idtft_hermitian"]


def fft_centered(x, fs=1):
    """
    Computes a *centered* FFT.

    Provides an FFT function where the output vector has the zero frequency
    at its center. This is more suitable for plotting.

    Parameters
    ----------
    x : array_like
        1-D sequence to compute the FFT upon

    Returns
    ------
    X : ndarray
        samples of the DFT of the input vector
    ff : ndarray
        vector of frequencies corresponding to the samples in X

    Other Parameters
    ----------------
    fs : real, optional
        sample frequency for the input vector (defaults to 1)
    """
    N = len(x)
    fs = float(fs)
    if np.mod(N, 2) == 0:
        k = np.arange(-N/2, N/2)  # N even
    else:
        k = np.arange(-(N-1)/2, (N-1)/2+1)  # N odd
    ff = k/(N/fs)
    X = sp.fftpack.fft(x)
    X = sp.fftpack.fftshift(X)
    return(ff, X)


def dtft(x, fs=1, t0=0):
    """
    Computes the discrete time Fourier transform (DTFT).

    Returns a function that is the DTFT of the given vector.

    Parameters
    ----------
    x :  array_like
        the 1-D vector to compute the DTFT upon

    Returns
    -------
    X : callable
        a function of frequency as in X(f), corresponding to the DTFT of x

    Other Parameters
    ----------------
    fs : real, optional
        sample frequency for the input vector (defaults to 1)
    t0 : real, optional
        the time when x[0] is sampled (defaults to 0). This is expressed
        in sample intervals.
    """
    return lambda f: np.sum(x*np.exp(-2j*np.pi*f/fs*(np.arange(len(x))-t0)))


def dtft_hermitian(x, fs=1):
    """
    Computes the discrete time Fourier transform of a hermitian vector.

    Since the input vector is hermitian, i.e. x[-n]=conjugate(x[n]), only
    a single side of it is sufficient for the computation of the DTFT,
    which is necessarily real.

    Parameters
    ----------
    x : array_like
        positive side of the 1-D sequence to compute the DTFT on.

    Returns
    -------
    X : callable
        a real function of frequency as in X(f), corresponding to the DTFT
        of x

    Other Parameters
    ----------------
    fs : real, optional
        sample frequency for the input vector (defaults to 1)
    """
    x = x.real
    return lambda f: \
        x[0]+np.sum(x[1:]*2*np.cos(2*np.pi*f/fs*np.arange(1, len(x))))


def _idtft(Ff, t, fs=1, **quad_opts):
    rr = sp.integrate.quad(lambda f: np.real(Ff(f*fs)*np.exp(2j*np.pi*f*t)),
                           -0.5, 0.5, **quad_opts)[0]
    ri = sp.integrate.quad(lambda f: np.imag(Ff(f*fs)*np.exp(2j*np.pi*f*t)),
                           -0.5, 0.5, **quad_opts)[0]
    return rr+ri


def idtft(Ff, tt, fs=1, **options):
    """Compute the inverse discrete time Fourier transform (IDTFT)

    Parameters
    ----------
    Ff : callable
        the function of frequency to compute the IDTFT from
    tt : real or array_like
        the time sample (or 1-D vector of) where to compute the IDTFT

    Returns
    -------
    x : ndarray
        the inverse discrete time Fourier transform at time values in tt

    Other Parameters
    ----------------
    fs : real, optional
        the sample frequency for the output sequence (defaults to 1)
    quad_xxx : various type
        Parameters prefixed by ``quad_`` are passed to the ``quad``
        function that is used internally as an integrator. Allowed options
        are ``quad_epsabs``, ``quad_epsrel``, ``quad_limit``, ``quad_points``.
        Do not use other options since they could break the integrator in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    See Also
    --------
    scipy.integrate.quad : integrator used internally.
        For the meaning of the integrator parameters.
    """
    # Manage optional parameters
    opts = idtft.default_options.copy()
    opts.update(options)
    quad_opts = {k[5:]: v for k, v in opts.iteritems()
                 if k.startswith('quad_')}
    # Do the computation
    if np.isscalar(tt):
        return _idtft(Ff, tt, fs, **quad_opts)
    else:
        return np.asarray([_idtft(Ff, t, fs, **quad_opts) for t in tt])

idtft.default_options = {'quad_epsabs': 1.49e-08,
                         'quad_epsrel': 1.49e-08,
                         'quad_limit': 50,
                         'quad_points': None}


def _idtft_hermitian(Ff, t, fs=1, **quad_opts):
    return 2*sp.integrate.quad(lambda f: np.real(Ff(f*fs))*np.cos(2*np.pi*f*t),
                               0, 0.5, **quad_opts)[0]


def idtft_hermitian(Ff, tt, fs=1, **options):
    """Compute the inverse discrete time Fourier transform (IDTFT) for a
    hermitian function of frequency.

    Since the input function is hermitian, i.e. F(-f)=conjugate(F(f)), it
    is sufficient to look at a single side of it. Furthermore, the IDTFT
    turns out to be a real vector.

    Parameters
    ----------
    Ff : callable
        the function of frequency to compute the IDTFT from. The values for
        a negative argument are not used.
    tt : real or array_like
        the time sample (or 1-D vector of) where to compute the IDTFT

    Returns
    -------
    x : ndarray
        the inverse discrete time Fourier transform at time values in tt

    Other Parameters
    ----------------
    fs : real, optional
        the sample frequency for the output sequence (defaults to 1)
    quad_xxx : various type
        Parameters prefixed by ``quad_`` are passed to the ``quad``
        function that is used internally as an integrator. Allowed options
        are ``quad_epsabs``, ``quad_epsrel``, ``quad_limit``, ``quad_points``.
        Do not use other options since they could break the integrator in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    See Also
    --------
    scipy.integrate.quad : integrator used internally.
        For the meaning of the integrator parameters.
    """
    # Manage optional parameters
    opts = idtft_hermitian.default_options.copy()
    opts.update(options)
    quad_opts = {k[5:]: v for k, v in opts.iteritems()
                 if k.startswith('quad_')}
    # Do the computation
    if np.isscalar(tt):
        return _idtft_hermitian(Ff, tt, fs, **quad_opts)
    else:
        return np.asarray([_idtft_hermitian(Ff, t, fs, **quad_opts)
                           for t in tt])

idtft_hermitian.default_options = idtft.default_options.copy()
