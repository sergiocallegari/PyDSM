# Copyright © 2013–2025 Sergio Callegari, Federico Bizzarri
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

"""Demo for LP ΔΣ modulator with low DC noise.

© 2013–2025 Sergio Callegari, Federico Bizzarri
All rights reserved.
"""


import numpy as np
from numpy.random import randn
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from pydsm.NTFdesign import ntf_fir_weighting
from pydsm.delsig import evalTF, dbv, dbp, simulateDSM

# Signal specification
fsig = 500.
B = 1000.
OSR = 64
fphi = B*OSR*2
# Lee constraint
H_inf = 1.5
# FIR Order for optimal NTF
order = 32


# Define the weighting function
# Start with a brickwall Lp function cutting at B
def w1(f):
    return 1. if f < B/fphi else 1e-12


# Introduce a weighting component that magnifies weight for
# frequencies < B/10
def weight2(f):
    return 50.-49.*f/(B/fphi/10.) if f < B/fphi/10. else 1.


# Get the final weighting function as the product of the two
def w2(f):
    return w1(f)*weight2(f)


print("... computing optimal NTF...")
ntf1 = ntf_fir_weighting(order, w1, quad_opts={'points': [B/fphi]})
ntf2 = ntf_fir_weighting(order, w2, quad_opts={'points': [B/fphi, B/fphi/10]})

# Prepare frequency axis for plotting
fmin = 10**np.ceil(np.log10(B/OSR/100))
fmax = fphi/2
ff = np.logspace(np.log10(fmin), np.log10(fmax), 1000)

resp_w1 = np.fromiter(map(w1, ff/fphi), np.double)
resp_ntf1 = np.abs(evalTF(ntf1, np.exp(1j*2*np.pi*ff/fphi)))
resp_w2 = np.fromiter(map(w2, ff/fphi), np.double)
resp_ntf2 = np.abs(evalTF(ntf2, np.exp(1j*2*np.pi*ff/fphi)))

ffa = np.logspace(np.log10(0.5E-5), np.log10(0.5), 1024)
vv_a2 = dbv(np.abs(evalTF(ntf2, np.exp(1j*2*np.pi*ffa))))

fig0 = plt.figure()
l = plt.plot(ff/fphi, dbp(resp_w1), 'b', label='on-off weighting')
plt.plot(ff/fphi, dbp(resp_w2), 'r', label='low-dc-noise weighting')
l = plt.plot(ff/fphi, dbv(resp_ntf1), 'b--', label='NTF - on-off')
plt.plot(ff/fphi, dbv(resp_ntf2), 'r--', label='NTF - low-dc-noise')
plt.xlim(1e-5, 1./2)
plt.xscale('log', base=10)
plt.gca().set_xticks([1E-5, 1./2])
plt.gca().set_xticklabels(['$10^{-5}$', r'$\frac{1}{2}$'])
plt.xlabel('$f$', x=1.)
plt.ylabel((r'$w(f)$, ' +
            r'$\left|\mathit{NTF}\,\left(\mathrm{e}^{\mathrm{i} 2\pi f}' +
            r'\right)\right|$ [dB]'), y=0.5)
plt.grid(True, 'both')
plt.suptitle('Noise weighting and NTF magnitude response')
plt.legend(loc='lower right')
plt.tight_layout(rect=[0, 0, 1, 0.98])

# Create an input signal for simulation
# Two signal components
fa = 53.
fb = 497.
periods = 2000
nfft_periods = 40
skip_periods = 10
tmax = int(np.ceil(periods*1./fa*fphi))
tplot0 = int(np.round(skip_periods*1./fa*fphi))
tplot2 = int(np.round((periods-2)*1./fa*fphi))
nn = np.arange(tmax+1)

# Three signal components available. Play freely with DC, slow and fast
# amplitudes
# Plot in paper obtained with
# c1=0 c2=0 c3=0.11111*np.sin(2*np.pi*fb*nn/fphi)*1E-3
#
# Constant component
c1 = 0.11111*0
# Slow component
c2 = 0*np.sin(2*np.pi*fa*nn/fphi)*1E-3
# Fast component
c3 = 0.11111*np.sin(2*np.pi*fb*nn/fphi)*1E-3
# dither
dither = 2.5E-6*randn(tmax+1)
sig = c1+c2+c3+dither

print("Simulating DSM")
ds1 = simulateDSM(sig, ntf1)[0]
ds2 = simulateDSM(sig, ntf2)[0]

ndt = tplot2-tplot0+1
dcval1 = np.sum(ds1[tplot0:tplot2])/ndt
dcval2 = np.sum(ds2[tplot0:tplot2])/ndt
print('Accuracy in the recovery of DC components')
print('recovered DC on-off', dcval1)
print('recovered DC low-dc-noise', dcval2)
print('original DC', np.sum(sig[tplot0:tplot2])/ndt)

print("Computing spectra")
NFFT = int(2**(np.ceil(np.log(nfft_periods*1./fa*fphi)/np.log(2))))
# This is a bit rough... there may be some residual power at the
# signal frequencies. To obtain a good plot make c1, c2, c3 very small
(psd1, freqs) = mlab.psd((ds1-sig)[tplot0:tplot2], Fs=2*fphi, NFFT=NFFT,
                         noverlap=NFFT//2)
(psd2, freqs) = mlab.psd((ds2-sig)[tplot0:tplot2], Fs=2*fphi, NFFT=NFFT,
                         noverlap=NFFT//2)
minfreq = 1
df = freqs[1]-freqs[0]
minidx = int(np.ceil(minfreq/df))

plt.figure()
plt.plot(freqs[minidx:], dbp(psd1[minidx:]), 'b', label='on-off')
plt.plot(freqs[minidx:], dbp(psd2[minidx:]), 'r', label='low-dc-noise')
plt.xlim(2, 128E3)
plt.xscale('log', base=10)
plt.xlabel('$f$', x=1.)
plt.ylabel('Quant. noise PSD [dB]')
plt.grid(True, 'both')
plt.legend(loc="upper left")
plt.tight_layout(rect=[0, 0, 1, 0.98])

plt.show()
