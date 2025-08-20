# Copyright © 2013–2025 Sergio Callegari
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

"""Demo code to reproduce the results shown in Figure 5 and in table 1,
of the paper:

S. Callegari 'Coding of Stereo Signals by a Single Digital ΔΣ Modulator'
IEEE 20th International Conference on Electronics, Circuits, and Systems
(ICECS 2013), pp. 589–592, Dec. 2013.

This code is about the proposed case, with a two band modulator and
multiplexing.

© 2013–2025 Sergio Callegari.
All rights reserved.
"""


import numpy as np
from scipy import signal
from pydsm.delsig import evalTF, dbv, dbp, simulateDSM, synthesizeNTF
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from dualNTF import mirroredNTF
from extract_tone import find_tone

# Set up problem data
# There are two signals, each of them with its own band.
# In this code, the reference case where an independent modulator is
#
B1 = 20E3     # Bandwidth of first signal
B2 = 20E3     # Bandwidth of second signal
OSR = 64      # Oversampling ratio
A1 = 0.2      # Amplitude of first signal
A2 = 0.44     # Amplitude of the second signal
f1 = 1003     # Frequency used to test the first test signal with a tone
f2 = 3193     # Frequency used to test the second test signal with a tone

Amax = 0.34   # Max amplitude used for testing

# Clock frequency of modulator
fphi = 2*OSR*(B1+B2)

# Cut off frequency of reconstruction filters
# ...set twice as wide as they should to better see artifacts
brk1 = B1/fphi
brk2 = 0.5-B2/fphi

# Create NTF with DELSIG synthesizeNTF
ntf_a = synthesizeNTF(order=4, osr=2*OSR, opt=3, H_inf=np.sqrt(1.5))
ntf_dual_a = mirroredNTF(ntf_a)

# Prepare frequency grid
ff_log = np.logspace(np.log10(0.5E-5), np.log10(0.5), 4096)
ff_lin = np.linspace(0, 0.5, 1024)

# Compute magnitude responses of NTF for plotting
ntf_a_mag = lambda f: np.abs(evalTF(ntf_a, np.exp(-2j*np.pi*f)))
vv_ntf_a_log = ntf_a_mag(ff_log)
vv_ntf_a_lin = ntf_a_mag(ff_lin)
ntf_dual_a_mag = lambda f: np.abs(evalTF(ntf_dual_a, np.exp(-2j*np.pi*f)))
vv_ntf_dual_a_log = ntf_dual_a_mag(ff_log)
vv_ntf_dual_a_lin = ntf_dual_a_mag(ff_lin)

# Design output filter
hz = signal.butter(4, 2*brk1, btype='low')

# Compute magnitude responses of output filter for plotting
hz_mag = lambda f: np.abs(evalTF(hz, np.exp(-2j*np.pi*f)))
vv_hz_log = hz_mag(ff_log)
vv_hz_lin = hz_mag(ff_lin)

print("Plotting the noise transfer function...")

# Plot magnitude responses with log scale - same as Fig. 5a
plt.figure()
plt.plot(ff_log, dbv(vv_ntf_dual_a_log))
plt.xlim(5E-5, 0.5)
plt.ylim(-119.99, 5)
plt.xscale('log', base=10)
plt.xlabel(r'$\hat f$', x=1.)
plt.ylabel(r'magnitude [dB]')
plt.gca().set_xticks([1E-5, 0.5])
plt.gca().set_xticklabels(['$10^{-5}$', r'$\frac{1}{2}$'])
plt.suptitle('NTF mag response - log scale')

# Plot magnitude responses with linear scale - same as Fig. 5b
plt.figure()
plt.plot(ff_lin, dbv(vv_ntf_dual_a_lin))
plt.xlim(0, 0.5)
plt.ylim(-100, 5)
plt.xlabel(r'$\hat f$', x=1.)
plt.ylabel(r'magnitude [dB]')
plt.gca().set_xticks([0, 1./2])
plt.gca().set_xticklabels(['$0$', r'$\frac{1}{2}$'])
plt.suptitle('NTF mag response - linear scale')

print("Doing time domain simulations...")

# Define test time for time domain simulations
Tstop = 200*fphi/1000  # approximately 200 periods at 1kHz
tt = np.arange(Tstop, dtype=int)
print("Test time:", len(tt), 'points')

# Generate signals for time domain simulations
print("Generate input signals ", end=" ")
s1 = A1*np.sin(2*np.pi*f1*tt/fphi)
print("* ", end="")
s2 = A2*np.sin(2*np.pi*f2*tt/fphi)
print("* ", end="")
clk = 2*(tt % 2)-1
print("# ", end="")
s2r = s2*clk
print("*")
dither_sigma = 1e-9  # Use a little dither
dither = np.random.randn(len(tt))*dither_sigma

# Do the simulation
print("Simulating DSM", end=" ")
xx1 = simulateDSM(s1+s2r+dither, ntf_dual_a)[0]
print("* ", end=""),
xx2 = xx1*clk
print("*")

# Plot a fragment of the modulator output - same as Fig. 5c
plt.figure()
plt.plot(tt[100:150], xx1[100:150], drawstyle='steps-post')
plt.xlabel(r'$n$', x=1.)
plt.ylabel(r'$x(nT)$', x=1.)
plt.ylim(-1.2, 1.2)
plt.suptitle('Fragment of modulator output - 1st signal')

# Reconstruct the signals
print("Applying a low pass filter")
yy1 = signal.lfilter(hz[0], hz[1], xx1)
yy2 = signal.lfilter(hz[0], hz[1], xx2)

# Plot fragment of first reconstructed signal - same as Figs. 5e, 5f
Tdisp_start = 10*fphi/1000
Tdisp_stop = 15*fphi/1000
ttp = np.arange(Tdisp_start, Tdisp_stop, dtype=int)
plt.figure()
plt.plot(ttp/fphi*1000, yy1[ttp], 'b', label='1st signal')
plt.plot(ttp/fphi*1000, yy2[ttp], 'r', label='2nd signal')
plt.xlabel(r'$t$ [ms]', x=1.)
plt.ylabel(r'$\hat u_1(t)$')
plt.suptitle('Fragment of reconstructed output')
plt.legend()

# Estimate some spectra
print("Compute some spectra (", end="")
Tsp_start = int(2*fphi/1000)
Tsp_stop = int(199*fphi/1000)
NFFT = 4096*64
(psd, freqs) = mlab.psd((xx1)[Tsp_start:Tsp_stop], Fs=fphi, NFFT=NFFT,
                        noverlap=NFFT//2, scale_by_freq=True)
print('@modulator, ', end="")
(psd1, freqs1) = mlab.psd((yy1)[Tsp_start:Tsp_stop], Fs=fphi, NFFT=NFFT,
                          noverlap=NFFT//2, scale_by_freq=True)
print('@1st signal, ', end="")
(psd2, freqs2) = mlab.psd((yy2)[Tsp_start:Tsp_stop], Fs=fphi, NFFT=NFFT,
                          noverlap=NFFT//2, scale_by_freq=True)
print('@2nd signal)')

print("Check overall power from psd")
print(" (should be 1): ", end="")
df = freqs[1]-freqs[0]
print(np.sum(psd)*df)

# Plot PDS of signal at modulator output
# This is the same as in Fig. 5b
# Try to express quantities in dBm (referred to 1 mW over a 50 Ohm resistor)
plt.figure()
plt.plot(freqs, dbp(psd)+30-10*np.log10(50))
plt.xlim(0, fphi/2)
plt.xlabel(r'$f$ [MHz]', x=1.)
plt.ylabel(r'PDS [dBm/Hz]')
plt.gca().set_xticks([0, fphi/2])
plt.gca().set_xticklabels(['$0$', r'$2.56$'])
plt.suptitle('PSD at modulator output')

# Plot PDS of first reconstructed signal
plt.figure()
plt.plot(freqs1, dbp(psd1))
plt.xscale('log', base=10)
plt.suptitle('PSD of reconstructed signal - 1st signal')

# Plot PDS of second reconstructed signal
plt.figure()
plt.plot(freqs2, dbp(psd2))
plt.xscale('log', base=10)
plt.suptitle('PSD of reconstructed signal - 2nd signal')

# Check modulator output when tone is subtracted
tone_xx1 = find_tone(xx1, f1, fphi, start=Tsp_start, N=Tsp_stop-Tsp_start,
                     output='polar')
tone_xx2 = find_tone(xx2, f2, fphi, start=Tsp_start, N=Tsp_stop-Tsp_start,
                     output='polar')
# Remove the tone, quick version
xx1q = xx1-tone_xx1[0]*np.cos(2*np.pi*f1/fphi*tt-tone_xx1[1])
xx2q = xx2-tone_xx2[0]*np.cos(2*np.pi*f2/fphi*tt-tone_xx2[1])
(psd1q, freqs1q) = mlab.psd(xx1q[Tsp_start:Tsp_stop], Fs=fphi, NFFT=NFFT,
                            noverlap=NFFT//2, scale_by_freq=True)
(psd2q, freqs2q) = mlab.psd(xx2q[Tsp_start:Tsp_stop], Fs=fphi, NFFT=NFFT,
                            noverlap=NFFT//2, scale_by_freq=True)

# Plot PDS of quantization noise for channel 1 at modulator output
# This is the same as Fig. 5g
plt.figure()
plt.plot(freqs1q/1000, dbp(psd1q)+30-10*np.log10(50))
plt.xscale('log', base=10)
plt.xlim(freqs1q[0]/1000, freqs1q[-1]/1000)
plt.xlabel(r'$f$ [kHz]', x=1.)
plt.ylabel(r'q. noise #1 [dBm/Hz]')
plt.suptitle('PSD of quantization noise - 1st signal')

# Plot PDS of quantization noise for channel 2 at modulator output
# This is the same as Fig. 5h
plt.figure()
plt.plot(freqs2q/1000, dbp(psd2q)+30-10*np.log10(50))
plt.xscale('log', base=10)
plt.xlim(freqs2q[0]/1000, freqs2q[-1]/1000)
plt.ylabel(r'q. noise #2 [dBm/Hz]')
plt.xlabel(r'$f$ [kHz]', x=1.)
plt.suptitle('PSD of quantization noise - 2nd signal')

# Compute merit factors referred to the modulator output
print('Compute merit factors for delta sigma sequences')
df1q = freqs1q[1]-freqs1q[0]
df2q = freqs2q[1]-freqs2q[0]
ifm1 = int(np.round((B1-freqs1q[0])/df1q)+1)
ifm2 = int(np.round((B2-freqs2q[0])/df2q)+1)
ibqnp1 = np.sum(psd1q[0:ifm1])*df1q
ibqnp2 = np.sum(psd2q[0:ifm1])*df2q
print('Noise floors', dbp(ibqnp1)+30-dbp(50), dbp(ibqnp2)+30-dbp(50), 'dBm')
print('SNRs', dbp(0.5*A1**2/ibqnp1), dbp(0.5*A2**2/ibqnp2), 'dB')
print('SNRmax', dbp(0.5*Amax**2/ibqnp1), dbp(0.5*Amax**2/ibqnp2), 'dB')

# Compute merit factors referred to the reconstructed signals
tone1 = find_tone(yy1, f1, fphi, start=Tsp_start, N=Tsp_stop-Tsp_start,
                  output='polar')
tone2 = find_tone(yy2, f2, fphi, start=Tsp_start, N=Tsp_stop-Tsp_start,
                  output='polar')
# Remove the tone
yy1q = yy1-tone1[0]*np.cos(2*np.pi*f1/fphi*tt-tone1[1])
yy2q = yy2-tone2[0]*np.cos(2*np.pi*f2/fphi*tt-tone2[1])
print('Signals at filter output:')
print('tone 1:', tone1)
print('tone 2:', tone2)

(psd1yq, freqs1yq) = mlab.psd(yy1q[Tsp_start:Tsp_stop], Fs=fphi, NFFT=NFFT,
                              noverlap=NFFT//2, scale_by_freq=True)
(psd2yq, freqs2yq) = mlab.psd(yy2q[Tsp_start:Tsp_stop], Fs=fphi, NFFT=NFFT,
                              noverlap=NFFT//2, scale_by_freq=True)
df1yq = freqs1yq[1]-freqs1yq[0]
df2yq = freqs2yq[1]-freqs2yq[0]
ibqnpy1 = np.sum(psd1yq)*df1yq
ibqnpy2 = np.sum(psd2yq)*df2yq
print('Noise floors', dbp(ibqnpy1)+30-dbp(50), dbp(ibqnpy2)+30-dbp(50), 'dB')
print('SNRs', dbp(0.5*tone1[0]**2/ibqnpy1), dbp(0.5*tone2[0]**2/ibqnpy2), 'dB')

# Plot PDS of quantization noise for channel 1 at filter output
# The same as in Fig. 5g, after filtering
plt.figure()
plt.plot(freqs1q/1000, dbp(psd1yq)+30-10*np.log10(50))
plt.xscale('log', base=10)
plt.xlim(freqs1q[0]/1000, freqs1q[-1]/1000)
plt.xlabel(r'$f$ [kHz]', x=1.)
plt.ylabel(r'q. noise #1 [dBm/Hz]')
plt.suptitle('PSD of quantization noise - 1st signal - filtered')

# Plot PDS of quantization noise for channel 2 at filter output
# The same as in Fig. 5h, after filtering
plt.figure()
plt.plot(freqs2q/1000, dbp(psd2yq)+30-10*np.log10(50))
plt.xscale('log', base=10)
plt.xlim(freqs2q[0]/1000, freqs2q[-1]/1000)
plt.ylabel(r'q. noise #2 [dBm/Hz]')
plt.xlabel(r'$f$ [kHz]', x=1.)
plt.suptitle('PSD of quantization noise - 2nd signal - filtered')

# Look for crosstalk
tonex1 = find_tone(yy2, f1, fphi, start=Tsp_start, N=Tsp_stop-Tsp_start,
                   output='polar')
tonex2 = find_tone(yy1, f2, fphi, start=Tsp_start, N=Tsp_stop-Tsp_start,
                   output='polar')
print('Crosstalk at filter output:')
print('tone 1:', tonex1, dbp(tonex1[0]**2)+30-dbp(50), 'dBm')
print('tone 2:', tonex2, dbp(tonex2[0]**2)+30-dbp(50), 'dBm')

plt.show()
