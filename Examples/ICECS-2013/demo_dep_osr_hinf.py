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

"""Demo code to replicate the result shown in Figure 4 of the paper:

S. Callegari 'Coding of Stereo Signals by a Single Digital ΔΣ Modulator'
IEEE 20th International Conference on Electronics, Circuits, and Systems
(ICECS 2013), pp. 589–592, Dec. 2013.

Copyright (c) Sergio Callegari, Federico Bizzarri 2012
All rights reserved.
"""

# Note: code rises a warning since in one of the cases being tested
# the NTF synthesis function is invoked with a too large maximum peak
# gain (Lee constraint) together with a very low modulator order.
# In this conditions, the synthesised NTF fails to reach the max gain.
# This does not affect the visualized trends, though.

import numpy as np
from pydsm.delsig import (synthesizeNTF, ds_f1f2, dbp)
from pydsm.NTFdesign.weighting import quantization_weighted_noise_gain
import matplotlib.pyplot as plt

# Set the OSRs and the h_inf (gamma) values used for the plots
# OSR values are arranged so that each one is twice the previous one
# hinf values are arranged so that each value is the square root of the
# previous one. Also set the orders used for testing.
osrs=np.asarray([16, 32, 64, 128, 256])
hinfs=np.asarray([2.25**(1/(2.**i)) for i in range(5)])
orders=np.asarray([1, 2, 3, 4, 5])

# Prepare a vector to store the noise gain for obtained at the many
# test conditions
g0s=np.zeros((orders.size, osrs.size))

# Synthesize the NTF for all the orders and OSR values at h_inf=1.5
for i, osr in enumerate(osrs):
    for j, order in enumerate(orders):
        ntf1=synthesizeNTF(order, osr=osr, opt=3, H_inf=1.5)
        f1, f2= ds_f1f2(osr)
        g0=quantization_weighted_noise_gain(ntf1, None, (f1, f2))
        # print order, osr, f1, f2, g0, dbp(g0)
        g0s[j, i]=g0
# Make the plot: quantization noise gain versus OSR for the different
# orders
plt.figure()
plt.xscale('log')
markers='ov^<>'
for i, order in enumerate(orders):
    plt.plot(osrs, dbp(g0s[i]), "-"+markers[i], label='%d' % order)
plt.legend(loc='upper right', fontsize=9)
plt.xlabel("OSR")
plt.ylabel("$P_N$ [dB]")
# Print some interesting data from the designs
print "Analysis of data for varying OSR"
for i, order in enumerate(orders):
    print "Order: %d. Average slope: %f dB per OSR doubling." % \
    (order, (dbp(g0s[i,-1])-dbp(g0s[i,0]))/(orders.size-1))

# Prepare a vector to store the noise gain for obtained at the many
# test conditions
g0s=np.zeros((orders.size, hinfs.size))

# Synthesize the NTF for all the orders and h_inf values at OSR=64
for i, hinf in enumerate(hinfs):
    for j, order in enumerate(orders):
        ntf1=synthesizeNTF(order, osr=64, opt=3, H_inf=hinf)
        f1, f2= ds_f1f2(osr)
        g0=quantization_weighted_noise_gain(ntf1, None, (f1, f2))
        #print order, hinf, f1, f2, g0, dbp(g0)
        g0s[j, i]=g0
# Make the plot: quantization noise gain versus h_inf for the different
# orders
plt.figure()
markers='ov^<>'
for i, order in enumerate(orders):
    plt.plot(np.arange(5)+1, dbp(g0s[i]), "-"+markers[i], label='%d' % order)
plt.legend(loc='lower right', fontsize=9)
plt.xlabel(r"$n$ so that $\gamma = \gamma_0^{\frac{1}{2^n}}$")
plt.ylabel("$P_N$ [dB]")
# Print some interesting data from the designs
print "\nAnalysis of data for varying h_inf"
for i, order in enumerate(orders):
    print "Order: %d. Average slope: %f dB per square root of h_inf." % \
    (order, (dbp(g0s[i,2])-dbp(g0s[i,1])))

plt.show()
