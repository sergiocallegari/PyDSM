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

"""Demo for ΔΣ modulator with brickwall LP filter.

© 2013–2025 Sergio Callegari, Federico Bizzarri
All rights reserved.
"""


import numpy as np
from pydsm.delsig import evalTF
from pydsm.delsig import dbv
import matplotlib.pyplot as plt
from pydsm.NTFdesign import ntf_fir_weighting
from pydsm.delsig import synthesizeNTF

# Signal specification
B = 1000.
OSR = 64
fphi = B*OSR*2
# Lee constraint
H_inf = 1.5
# FIR Order for optimal NTF
order = 10
# Order for reference NTF (DELSIG)
delsig_order = 2


# On-off weighting function
def w1(f):
    return 1. if f < B/fphi else 1e-12


print("... plotting filter magnitude response")
ff = np.logspace(np.log10(0.5E-5), np.log10(0.5), 1024).reshape((-1, 1))
vv = np.fromiter(map(w1, ff), np.double)
plt.figure()
plt.plot(ff, vv)
plt.xlim(1e-5, 1./2)
plt.ylim(0, 1.1)
plt.xscale('log', base=10)
plt.gca().set_xticks([1E-5, 1./2])
plt.gca().set_xticklabels(['$10^{-5}$', r'$\frac{1}{2}$'])
plt.xlabel('$f$', x=1.)
plt.ylabel('$w(f)$', y=0.9)
plt.grid(True, 'both')
plt.suptitle('Weighting function')
plt.tight_layout(rect=[0, 0, 1, 0.98])

print("... computing optimal NTF")
ntf_opti = ntf_fir_weighting(order, w1, H_inf=H_inf)
ntf_opti_mag = lambda f: np.abs(evalTF(ntf_opti, np.exp(-2j*np.pi*f)))

print("... computing delsig NTF")
ntf_delsig = synthesizeNTF(delsig_order, OSR, 3, H_inf, 0)
ntf_delsig_mag = lambda f: np.abs(evalTF(ntf_delsig, np.exp(-2j*np.pi*f)))

print("... plotting optimal NTF amplitude response")
vv_mag = dbv(ntf_opti_mag(ff))
vv_mag_delsig = dbv(ntf_delsig_mag(ff))
plt.figure()
plt.plot(ff, vv_mag, label='Proposed')
plt.plot(ff, vv_mag_delsig, 'r-o', linewidth=0.5, markevery=24, markersize=3,
         label='Reference')
plt.xlim(1e-5, 1./2)
# plt.ylim(0,1.1)
plt.xscale('log', base=10)
plt.gca().set_xticks([1E-5, 1./2])
plt.gca().set_xticklabels(['$10^{-5}$', r'$\frac{1}{2}$'])
plt.xlabel('$f$', x=1.)
plt.ylabel((r'$\left|\mathit{NTF}\,\left(\mathrm{e}^{\mathrm{i} 2\pi f}' +
            r'\right)\right|$ [dB]'), y=0.59)
plt.grid(True, 'both')
plt.suptitle('Noise transfer function (magnitude)')
plt.legend(loc="lower right")
plt.tight_layout()

plt.show()
