# -*- coding: utf-8 -*-

# Copyright (c) 2013, Sergio Callegari, Federico Bizzarri
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

u"""Demo for LP DS modulator with psychoacoustic noise weighting.

Copyright © 2013 Sergio Callegari, Federico Bizzarri
All rights reserved.
"""

from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from pydsm.delsig import dbv, dbp
from pydsm.NTFdesign import ntf_dunn, ntf_fir_audio_weighting
from pydsm.delsig import simulateDSM, evalTF
import matplotlib.mlab as mlab
from pydsm.audio_weightings import f_weighting

# Signal specification
fsig = 1000.
B = 20500.
osr = 64
fphi = B*osr*2
# Lee constraint
H_inf = 1.5
# FIR Order
order = 30

# Computing the optimal NTF and the reference NTF
print("... computing optimal NTF")
# Reducing the feastol may help the optimization converge.
# Until reltol is very strict, subobtima are often returned
# Using a large normalization constant may help getting closer to the real
# optimum but may also hinder convergence
opti_ntf = ntf_fir_audio_weighting(order, osr, f_weighting,
                                   H_inf=H_inf,
                                   normalize=1E3,
                                   cvxopt_opts={'reltol': 1E-12,
                                                'abstol': 1E-10,
                                                'feastol': 1E-2})

dunn_ntf = ntf_dunn(3, osr, H_inf)

fmin = 10
fmax = fphi/2
ff = np.logspace(np.log10(fmin), np.log10(fmax), 1000)

resp_w = f_weighting(ff)
resp_opti = np.abs(evalTF(opti_ntf, np.exp(1j*2*np.pi*ff/fphi)))
resp_ref = np.abs(evalTF(dunn_ntf, np.exp(1j*2*np.pi*ff/fphi)))

# First figure. This provides the weighting function and the NTFs
fig0 = plt.figure()
plt.plot(ff/fphi, dbp(resp_w), 'b', label='audio weighting')
plt.plot(ff/fphi, dbv(resp_opti), 'r', label='proposed NTF')
l = plt.plot(ff/fphi, dbv(resp_ref), 'g', label="Dunn's NTF")
plt.suptitle('Audio weighting and NTF magnitude response')
plt.xlim(10E-5, 1./2)
plt.ylim(-140, 20)
plt.xscale('log', basex=10)
plt.gca().set_xticks([1E-5, 1./2])
plt.gca().set_xticklabels(['$10^{-5}$', r'$\frac{1}{2}$'])
plt.xlabel('$f$', x=1.)
plt.ylabel((r'$w(f)$, ' +
            r'$\left|\mathit{NTF}\,\left(\mathrm{e}^{\mathrm{i} 2\pi f}' +
            r'\right)\right|$ [dB]'), y=0.5)
plt.grid(True, 'both')
plt.legend(loc='lower right')
plt.tight_layout(rect=[0, 0, 1, 0.98])

# Start and stop time for DS simulation
Tstop = 400000
Tstart = 40000
dither_sigma = 0
# Use very little signals, otherwise it is troublesome to
# remove them completely to evaluate the quantization noise
A = 0.01
DT = 4*fphi/fsig

# Setting up DS simulation
tt = np.asarray(range(int(Tstop)))
uu = A*np.sin(2*np.pi*fsig/fphi*tt)
dither = np.random.randn(len(uu))*dither_sigma
uud = uu+dither

# Simulating the DS
# Some of the commented bits are used to obtain the quantizazion noise pdf
print("Simulating optimal NTF")
xx_opti = simulateDSM(uud, opti_ntf)[0]
xx_ref = simulateDSM(uud, dunn_ntf)[0]

NFFT = 4096*32
# Possibly this can be improved by extracting the tone from the
# output to do the subtraction
(psd1, freqs) = mlab.psd((xx_opti-uu)[Tstart:Tstop], Fs=fphi, NFFT=NFFT,
                         noverlap=NFFT/2, scale_by_freq=True)
(psd2, reqs) = mlab.psd((xx_ref-uu)[Tstart:Tstop], Fs=fphi, NFFT=NFFT,
                        noverlap=NFFT/2, scale_by_freq=True)

minfreq = 1
df = freqs[1]-freqs[0]
print('frequency step in psd data=', df)

psd1 *= df
psd2 *= df

minidx = int(np.floor(minfreq/df))
print('minimum frequency index for plots', minidx)
maxidx = int(np.ceil(22E3/df))

# Here we plot the spectra obtained by time domain data
plt.figure()
plt.semilogx(freqs[minidx:], dbp(psd1[minidx:]), 'r', label='proposed')
plt.semilogx(freqs[minidx:], dbp(psd2[minidx:]), 'g', label="Dunn's",
             alpha=0.4)
plt.semilogx(ff, dbp(1.3*resp_opti**2/(fphi)), 'k', label='proposed NTF')
plt.legend(loc='lower right')
plt.title('Spectra from time domain data')
plt.tight_layout(rect=[0, 0, 1, 0.98])

ww = np.zeros_like(psd1)
ww = np.asarray([max(f_weighting(freq), 1E-200) for freq in freqs])

plt.figure()
plt.plot(freqs[minidx:maxidx], dbp((psd1*ww)[minidx:maxidx]), 'r',
         label='proposed')
l = plt.plot(freqs[minidx:maxidx], dbp((psd2*ww)[minidx:maxidx]), 'g',
             label="Dunn's")
l[0].set_dashes([2, 2])
plt.xscale('log', basex=10)
plt.xlim(freqs[minidx], 22E3)
plt.ylim(-240, -120)
plt.xlabel('$f$', x=1.)
plt.ylabel('Loudness dB', y=0.5)
plt.grid(True, 'both')
plt.title("Perceived noise level")
plt.legend(loc='lower center')
plt.tight_layout(rect=[0, 0, 1, 0.98])

plt.show()
