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

u"""Demo for multiband DS modulator with brickwall filter.

Copyright © 2013 Sergio Callegari, Federico Bizzarri
All rights reserved.
"""

import numpy as np
from pydsm.delsig import evalTF
from pydsm.delsig import dbv
import matplotlib.pyplot as plt
from pydsm.NTFdesign.weighting import synthesize_ntf_from_noise_weighting

import matplotlib
matplotlib.rc('font', size=16, family='serif')

# Signal specification
B=1000.
OSR=64
fphi=B*OSR*2
# Lee constraint
H_inf=1.5
# FIR Order for optimal NTF
order=32

# On-off weighting function
def w1(f):
    return 1. if f < B/fphi or 8*B/fphi < f < 15*B/fphi else 1e-12


print("... plotting filter magnitude response")
ff=np.logspace(np.log10(0.5E-5),np.log10(0.5),1024).reshape((-1,1))
vv=np.asarray(map(w1,ff))
plt.figure()
plt.plot(ff,vv)
plt.xlim(1e-5,1./2)
plt.ylim(0,1.1)
plt.xscale('log',basex=10)
plt.gca().set_xticks([1E-5,1./2])
plt.gca().set_xticklabels(['$10^{-5}$',r'$\frac{1}{2}$'])
plt.xlabel('$f$', x=1.)
plt.ylabel('$w(f)$', y=0.9)
plt.grid(True,'both')
plt.suptitle('Weighting function')
plt.tight_layout(rect=[0,0,1,0.98])

print("... computing optimal NTF")
ntf_opti=synthesize_ntf_from_noise_weighting(order, w1, H_inf=H_inf,
                                             options={'reltol':1E-11,
                                                      'abstol':1E-9,
                                                      'feastol':1E-9})
ntf_opti_mag=lambda f: np.abs(evalTF(ntf_opti,np.exp(-2j*np.pi*f)))


print("... plotting optimal NTF amplitude response")
vv_mag=dbv(ntf_opti_mag(ff))
plt.figure()
plt.plot(ff,vv_mag)
plt.xlim(1e-5,1./2)
#plt.ylim(0,1.1)
plt.xscale('log',basex=10)
plt.gca().set_xticks([1E-5,1./2])
plt.gca().set_xticklabels(['$10^{-5}$',r'$\frac{1}{2}$'])
plt.xlabel('$f$', x=1.)
plt.ylabel((r'$\left|\mathit{NTF}\,\left(\mathrm{e}^{\mathrm{i} 2\pi f}' +
            r'\right)\right|$ [dB]'), y=0.59)
plt.grid(True,'both')
plt.suptitle('Noise transfer function (magnitude)')
plt.tight_layout()


plt.show()