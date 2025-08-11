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

u"""
Psychoacoustic NTF design (:mod:`pydsm.NTFdesign.psychoacoustic`)
=================================================================

This modules provides a psychoacoustic aware design method for the NTF
of ΔΣ modulators used in audio applications.

.. currentmodule:: pydsm.NTFdesign.psychoacoustic


Functions
---------

.. autosummary::
   :toctree: generated/

    dunn_optzeros            -- helper function for normalized NTF zeros
    dunn_optzeros_cplx       -- helper function for NTF zeros
    ntf_dunn                 -- NTF design a la Dunn
    ntf_fir_audio_weighting  -- FIR NTF design with psychoacustic weighting


Deprecated functions
--------------------

.. autosummary::

   synthesize_ntf_dunn                 -- alias of ntf_dunn
   synthesize_ntf_from_audio_weighting -- alias of ntf_fir_audio_weighting
"""

from __future__ import division, print_function

import numpy as np
from .delsig import synthesizeNTF
from ..delsig import undbp as undbp
from .. import audio_weightings
from .weighting import ntf_fir_weighting
from warnings import warn
from ..exceptions import PyDsmDeprecationWarning

__all__ = ["dunn_optzeros", "dunn_optzeros_cplx", "synthesize_ntf_dunn",
           "synthesize_ntf_from_audio_weighting",
           "ntf_dunn", "ntf_fir_audio_weighting"]


def dunn_optzeros(n):
    """
    Helper function for synthesising psychoacoustically optimal modulators.

    Returns the (normalized) zeros which minimize the in-band noise power of
    a delta-sigma modulator's NTF after weighting by the F-weighting.

    Parameters
    ----------
    n : int
        the number of optimized zeros to return

    Returns
    -------
    zeros : ndarray of reals
        the zeros for the modulator as normalized angular frequencies.

    Notes
    -----
    The zeros are always located on the complex unit circle. In other words,
    the zeros are returned as frequencies. For homogeneity with DELSIG's
    ds_optzeros, the frequencies are normalized with respect to the signal
    bandwidth, that is fixed at 22.05 kHz.

    This function is the equivalent of ds_optzeros in DELSIG. The tabled
    zeros delivered by this function are from [1]_.  Note that this function
    does not return the values that are tabled in [1]_, but scales them
    by the reference audio bandwidth used in [1]_, namely 22.05 kHz.

    .. [1] Chris Dunn and Mark Sandler, "Psychoacoustically Optimal
       Sigma Delta Modulation," J. Audio Eng. Soc., Vol. 45, No. 4, pp.
       212 - 223 (1997 April)
    """
    # These are the optimal zero placements in kHz for a 22.05 kHz bandwidth
    # as found in the paper by Dunn and Sandler
    zero_freqs_unnorm = [[0.0],
                         [4.014, -4.014],
                         [0.0, 6.443, -6.443],
                         [3.590, -3.590, 11.954, -11.954],
                         [0.0, 4.308, -4.308, 12.959, -12.959],
                         [3.325, -3.325, 7.078, -7.078, 13.389, -13.389],
                         [0.0, 0.0, 4.017, -4.017, 10.471, -10.471, 13.842,
                          -13.842],
                         [2.933, -2.933, 5.167, -5.167, 12.012, -12.012,
                          14.381, -14.381]]
    if n > 8:
        raise ValueError('Optimized zeros for n>14 are not available.')
    return np.asarray(zero_freqs_unnorm[n-1])/22.05


def dunn_optzeros_cplx(order, osr):
    """
    Helper function for synthesising psychoacoustically optimal modulators.

    Returns the complex zeros which minimize the in-band noise power of
    a delta-sigma modulator's NTF after weighting by the F-weighting.
    The signal bandwidth is 22.05 kHz.

    Parameters
    ----------
    n : int
        the number of optimized zeros to return
    osr : real
        the oversampling ratio

    Returns
    -------
    zeros : ndarray of complex
        the zeros for the modulator as complex values

    See Also
    --------
    dunn_optzeros : Dunn's optimal zeros
    """
    w2 = dunn_optzeros(order)/osr*np.pi
    return np.exp(1j*w2)


def ntf_dunn(order=3, osr=64, H_inf=1.5):
    """
    Synthesizes an NTF for a DS audio modulator by Dunn's approach.

    The signal bandwidth is 22.05 kHz.

    Parameters
    ----------
    order : int, optional
        the order of the modulator, defaults to 3
    osr : float, optional
        the oversamping ratio (based on the actual signal bandwidth)
    H_inf : real, optional
        max allowed peak value of the NTF. Defaults to 1.5

    Returns
    -------
    ntf : tuple
        noise transfer function in zpk form.

    Warns
    -----
    As in DELSIG's synthesizeNTF

    Notes
    -----
    This is not exactly Dunn's method (but it should be equivalent or
    slightly better). In fact, to avoid re-implementing the pole selection
    algorithm, that in [1]_ is based on a Butterworth synthesis, here
    the pole selection logic used in DELSIG's synthesizeNTF is recycled.
    This should not make a big difference since DELSIG logic is anyway based
    on a pole positioning aimed at obtaining a maximally flat response of the
    NTF denominator in the signal band. This has the advantage of
    automatically controlling the peak gain of the NTF, as in the Lee
    criterion.

    Parameter H_inf is used to enforce the Lee stability criterion.

    .. [1] Chris Dunn and Mark Sandler, "Psychoacoustically Optimal
       Sigma Delta Modulation," J. Audio Eng. Soc., Vol. 45, No. 4, pp.
       212 - 223 (1997 April)

    See Also
    --------
    pydsm.delsig.synthesizeNTF : DELSIG's optimal NTF design strategy.
    """
    return synthesizeNTF(order, osr, dunn_optzeros_cplx(order, osr), H_inf, 0)


def ntf_fir_audio_weighting(
        order, osr,
        audio_weighting='f',
        audio_band=22.05E3,
        max_attn=120,
        H_inf=1.5,
        normalize="auto", **options):
    u"""
    Synthesize a FIR NTF based on an audio weighting function.

    The ΔΣ modulator NTF is designed after an audio weigthing function stating
    how loudly noise is perceived at the various frequencies.

    Parameters
    ----------
    order : int
        Delta sigma modulator order
    osr : float
        the oversampling ratio
    audio_weighting : callable or string
        audio weighting function. This is a function taking a frequency and
        expressing the weighting at that frequency in terms of acoustic power.
        Functions in the :mod:`pydsm.audio_weightings` module are suitable
        here. Note that the function argument frequency is a real frequency in
        Hz.
    audio_band : float, optional
        how large the audio bandwidth to consider. The signal band is from
        0 to audio_band Hz. Defaults to 22.05 kHz
    max_attn : float, optional
        clip very large attenuations to this value (in dB). This helps the
        convergenze of the optimization routine used to design the NTF.
        Defaults to 120 dB.
    H_inf : real, optional
        Max peak NTF gain, defaults to 1.5, used to enforce the Lee criterion
    normalize : string or real, optional
        Normalization to apply to the quadratic form used in the NTF
        selection. Defaults to 'auto' which means setting the top left entry
        in the matrix Q defining the quadratic form to 1.

    Returns
    -------
    ntf : ndarray
        FIR NTF in zpk form

    Other parameters
    ----------------
    show_progress : bool, optional
        provide extended output, default is True and can be updated by
        changing the function ``default_options`` attribute.
    cvxopt_opts : dictionary, optional
        A dictionary of options for the ``cvxopt`` optimizer
        Allowed options include:

        ``maxiters``
            Maximum number of iterations (defaults to 100)
        ``abstol``
            Absolute accuracy (defaults to 1e-7)
        ``reltol``
            Relative accuracy (defaults to 1e-6)
        ``feastol``
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break :mod:`cvxopt` in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.
    quad_opts : dictionary, optional
        Parameters to be passed to the :func:`scipy.minimize.quad` function
        used internally as an integrator. Allowed options are ``epsabs``,
        ``epsrel``, ``limit``, ``points``. Do not use other options since they
        could break the integrator in unexpected ways. Defaults can be set by
        changing the function ``default_options`` attribute.

    See Also
    --------
    scipy.integrate.quad : for the meaning of the integrator parameters.
    cvxopt : for the optimizer parameters
    """
    if not callable(audio_weighting):
        if audio_weighting == 'f':
            audio_weighting = audio_weightings.f_weighting
        elif audio_weighting == 'a':
            audio_weighting = audio_weightings.a_weighting
        elif audio_weighting == 'b':
            audio_weighting = audio_weightings.b_weighting
        elif audio_weighting == 'c':
            audio_weighting = audio_weightings.c_weighting
        elif audio_weighting == 'd':
            audio_weighting = audio_weightings.d_weighting

    def w(f):
        ma = undbp(-max_attn)
        fx = f*audio_band*2*osr
        w = audio_weighting(fx) if fx <= audio_band else 0
        return max(w, ma)

    return ntf_fir_weighting(order, w, H_inf, normalize, **options)

ntf_fir_audio_weighting.default_options = \
    ntf_fir_weighting.default_options.copy()


# Following part is deprecated

def synthesize_ntf_dunn(order=3, osr=64, H_inf=1.5):
    """
    Alias of :func:`ntf_dunn`

    .. deprecated:: 0.11.0
       Function has been moved to the :mod:`NTFdesign` module with
       name :func:`ntf_dunn`.
    """
    warn("Function superseded by ntf_dunn in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_dunn(order, osr, H_inf)


def synthesize_ntf_from_audio_weighting(
        order, osr,
        audio_weighting=audio_weightings.f_weighting,
        audio_band=22.05E3,
        max_attn=120,
        H_inf=1.5,
        normalize="auto", **options):
    """
    Alias of :func:`ntf_fir_audio_weighting`

    .. deprecated:: 0.11.0
        Function has been moved to the :mod:`NTFdesign` module with
        name :func:`ntf_fir_audio_weighting`.
    """
    warn("Function superseded by ntf_fir_audio_weighting in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_fir_audio_weighting(order, osr, audio_weighting, audio_band,
                                   max_attn, H_inf, normalize, **options)
