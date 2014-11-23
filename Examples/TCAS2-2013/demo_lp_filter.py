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

u"""Demo for DS modulator with assigned reconstruction filter.

Copyright © 2013 Sergio Callegari, Federico Bizzarri
All rights reserved.
"""

from __future__ import division, print_function

import numpy as np
from scipy.signal import butter
from pydsm.delsig import evalTF, synthesizeNTF
from pydsm.delsig import dbv, dbp
from pydsm.NTFdesign.filter_based import quantization_noise_gain
import matplotlib.pyplot as plt
from pydsm.NTFdesign.filter_based import synthesize_ntf_from_filter

# Signal specification
B = 1000.
OSR = 64
fphi = B*OSR*2
# Lee constraint
H_inf = 1.5
# FIR Order for optimal NTF
order = 10
# Order for reference NTF (DELSIG)
delsig_order = 4

# Define the reconstruction filter
# Cuts at 2B in order to have negligible attenuation in the signal band
hz = butter(1, B/fphi, btype='low', output='zpk')
hz_mag = lambda f: np.abs(evalTF(hz, np.exp(-2j*np.pi*f)))


print("... plotting filter magnitude response")
ff = np.logspace(np.log10(0.5E-5), np.log10(0.5), 1024).reshape((-1, 1))
vv = dbv(np.abs(hz_mag(ff)))
plt.figure()
plt.plot(ff, vv)
plt.xlim(1e-5, 1./2)
plt.ylim(-50, 5)
plt.xscale('log', basex=10)
plt.gca().set_xticks([1E-5, 1./2])
plt.gca().set_xticklabels(['$10^{-5}$', r'$\frac{1}{2}$'])
plt.xlabel('$f$', x=1.)
plt.ylabel('$w(f)$ [dB]', y=0.9)
plt.grid(True, 'both')
plt.suptitle('Output filter magnitude response')
plt.tight_layout(rect=[0, 0, 1, 0.98])

print("... computing optimal NTF")
ntf_opti = synthesize_ntf_from_filter(order, hz, H_inf=H_inf)
ntf_opti_mag = lambda f: np.abs(evalTF(ntf_opti, np.exp(-2j*np.pi*f)))

print("... computing reference NTF")
ntf_delsig = synthesizeNTF(delsig_order, OSR, 3, H_inf, 0)
ntf_delsig_mag = lambda f: np.abs(evalTF(ntf_delsig, np.exp(-2j*np.pi*f)))

print("... plotting optimal NTF amplitude response")
vv_mag = dbv(ntf_opti_mag(ff))
vv_ref = dbv(ntf_delsig_mag(ff))
plt.figure()
plt.plot(ff, vv_mag, 'b', label='proposed')
l = plt.plot(ff, vv_ref, 'r', label='reference')
plt.xlim(1e-5, 1./2)
plt.ylim(-100, 5)
plt.xscale('log', basex=10)
plt.gca().set_xticks([1E-5, 1./2])
plt.gca().set_xticklabels(['$10^{-5}$', r'$\frac{1}{2}$'])
plt.xlabel('$f$', x=1.)
plt.ylabel((r'$\left|\mathit{NTF}\,\left(\mathrm{e}^{\mathrm{i} 2\pi f}' +
            r'\right)\right|$ [dB]'), y=0.59)
plt.grid(True, 'both')
plt.suptitle('Noise transfer function (magnitude)')
plt.legend(loc="lower right")
plt.tight_layout(rect=[0, 0, 1, 0.98])

print("... computing some interesting data")
p1 = quantization_noise_gain(ntf_opti, hz)
p2 = quantization_noise_gain(ntf_delsig, hz)

print("Overall noise power amplification/attenuation")
print("Proposed", p1, "(", dbp(p1), "dB)")
print("Reference", p2, "(", dbp(p2), "dB)")

plt.show()
