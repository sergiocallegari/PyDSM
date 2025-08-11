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


import numpy as np
from scipy import signal

__all__ = ['find_tone']


def find_tone(x_in, f, fs=1, start=0, N=None, adjust_N=False, window='hamming',
              output='complex'):
    """Extracts a tone out of a signal where it might be embedded in noise.

    This is a complicate function to extract information (amplitude and
    phase) about a tone embedded in noise. In spite of its sophistication
    this function may not always work correctly.

    Parameters
    ----------
    x_in : array of floats
        input signal
    f : float
        frequency of tone being looked for
    fs : float
        sample frequency in the input vector x_in.
        Defaults to 1.
    start : integer
        first input sample to consider.
        Defaults to 0.
    N : integer or None
        number of input samples to consider. If it is None, then the number
        of samples taken in consideration is autoadjusted, using all the
        available samples, starting from start.
        Defaults to None.
    adjust_N : bool
        whether to adjust the number of samples to consider, in order to
        get an integer number of periodos of the tone being searched.
        Defaults to False.
    window : string, float or tuple
        type of window function used before practicing a Fourier Transform
        to look for the tone. Can be any window supported by the function
        get_window in the scipy signal module.
        Defaults to 'hanning'.
    output : string
        type of output to provide. Can be 'complex' or 'xy'.
        Defaults to 'complex'.

    Returns
    -------
    tone : complex or couple of floats.
        output format is controlled by input parameter `output`.
        If format is 'complex', then ouput encodes amplitude and phase of the
        input tone as a complex exponential. If format is 'xy' then real and
        'imaginary part' are provided separately.
    """
    M = len(x_in)
    fn = float(f)/fs
    samples_per_period = 1/fn
    Nmax = M-start
    if N is None:
        N = Nmax
    Pmax = int(np.floor(Nmax/samples_per_period))
    P = np.min((Pmax, int(np.round(N/samples_per_period))))
    if adjust_N:
        N2 = int(np.round(P*samples_per_period))
    else:
        N2 = int(np.min((Nmax, N)))
    if window == 'boxcar' or window is None:
        data = x_in[start:N2+start]
    else:
        win = signal.get_window(window, N2)
        adj = np.sum(win)/N2
        win = win/adj
        data = x_in[start:N2+start]*win
    x_cos = np.cos(2*np.pi*fn*np.arange(start, N2+start))
    x_sin = np.sin(2*np.pi*fn*np.arange(start, N2+start))
    b = 2*np.dot(x_cos, data)/N2
    a = 2*np.dot(x_sin, data)/N2
    if output == 'complex':
        return b+1j*a
    elif output == 'xy':
        return (b, a)
    else:
        return (np.sqrt(a**2+b**2), np.arctan2(a, b))
