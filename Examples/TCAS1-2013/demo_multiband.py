# -*- coding: utf-8 -*-

# Copyright (c) Sergio Callegari, Federico Bizzarri 2012
# All rights reserved.

"""Demo for DS modulator with multiband filter.

Copyright (c) Sergio Callegari, Federico Bizzarri 2012
All rights reserved.
"""

import numpy as np
import scipy as sp
__import__("scipy.signal")
import matplotlib.pyplot as plt
from pydsm.ir import impulse_response
from pydsm.delsig import simulateDSM, evalTF
from pydsm.delsig import dbv, dbp
from pydsm.NTFdesign.filter_based import quantization_noise_gain, \
    synthesize_ntf_from_filter_ir

# Signal specification
fsig1=1000.
fsig2=10000.
B1=400.
B2=4000.
OSR=64
fphi=(B1+B2)*OSR*2
# Lee constraint
H_inf=1.5
# FIR Order
order=50
# Signal amplitude
A1=0.45
A2=0.45

# Generate filter. Transfer function is normalized to be 0dB in pass band
print("...generating filter")
# Care: in butter the cut of frequency is specified as a number from 0 to 1
# where 1 is fphi/2, not fphi
w01=2*fsig1/fphi
B01=2*B1/fphi
w11=(np.sqrt(B01**2+4*w01**2)-B01)/2
w21=(np.sqrt(B01**2+4*w01**2)+B01)/2
hz1=sp.signal.butter(2, [w11,w21], 'bandpass', output='zpk')
w02=2*fsig2/fphi
B02=2*B2/fphi
w12=(np.sqrt(B02**2+4*w02**2)-B02)/2
w22=(np.sqrt(B02**2+4*w02**2)+B02)/2
hz2=sp.signal.butter(2, [w12,w22], 'bandpass', output='zpk')
hz1_ab=sp.signal.zpk2tf(*hz1)
hz2_ab=sp.signal.zpk2tf(*hz2)
hz_ab_d=np.polymul(hz1_ab[1],hz2_ab[1])
hz_ab_n1=np.polymul(hz1_ab[0],hz2_ab[1])
hz_ab_n2=np.polymul(hz2_ab[0],hz1_ab[1])
hz_ab_n=np.polyadd(hz_ab_n1,hz_ab_n2)
hz=sp.signal.tf2zpk(hz_ab_n,hz_ab_d)

# Compute impulse response
print("...computing impulse response of filter")
hz_ir=impulse_response(hz, db=60)

# Compute the optimal NTF
print("... computing optimal NTF")
ntf_opti=synthesize_ntf_from_filter_ir(order, hz_ir, H_inf=H_inf)
ntf_opti_zpk=(np.roots(ntf_opti),np.zeros(order),1)

# Determine freq values for which plots are created
fmin=10**np.ceil(np.log10(10))
fmax=10**np.floor(np.log10(fphi/2))
ff=np.logspace(np.log10(fmin),np.log10(fmax),1000)

# Compute frequency response data
resp_filt=np.abs(evalTF(hz,np.exp(1j*2*np.pi*ff/fphi)))
resp_opti=np.abs(evalTF((ntf_opti,[1]),np.exp(1j*2*np.pi*ff/fphi)))

# Plot frequency response
plt.figure()
plt.semilogx(ff,dbv(resp_filt), 'b', label="Output filter")
plt.semilogx(ff,dbv(resp_opti), 'r', label="Optimal NTF")
plt.legend(loc="lower right")
plt.suptitle("Output filter and NTFs")

# Check merit factors
ffl=np.linspace(fmin, fmax, 1000)
pg_opti=np.abs(evalTF((ntf_opti,[1]),np.exp(1j*2*np.pi*ffl/fphi)))* \
    np.abs(evalTF(hz,np.exp(1j*2*np.pi*ffl/fphi)))
plt.figure()
plt.plot(ffl,pg_opti**2,'r', label="Optimal NTF")
#plt.legend(loc="upper right")
plt.suptitle("Merit factor integrand")

# Compute expected behavior
sigma2_e=1./3
noise_power_opti_1=quantization_noise_gain(hz, (ntf_opti,[1]))*sigma2_e
print("Expected optimal noise level {} ({} dB)\nExpected SNR {} dB".format( \
    noise_power_opti_1, dbp(noise_power_opti_1), \
        dbp(0.5*A1**2+0.5*A2**2)-dbp(noise_power_opti_1)))

# Start and stop time for DS simulation
Tstop=100E3
Tstart=40E3
dither_sigma=1e-6

# Set up DSM simulation
tt=np.asarray(xrange(int(Tstop)))
uu=A1*np.sin(2*np.pi*fsig1/fphi*tt)+A2*np.sin(2*np.pi*fsig2/fphi*tt)
dither=np.random.randn(len(uu))*dither_sigma
uud=uu+dither

# Simulate the DSM
print("Simulating optimal NTF")
xx_opti = simulateDSM(uud, ntf_opti_zpk)[0]

print("Applying the reconstrution filter")
hz_ab=sp.signal.zpk2tf(*hz)
uu_filt=sp.signal.lfilter(hz_ab[0],hz_ab[1],uu)
xx_opti_filt=sp.signal.lfilter(hz_ab[0],hz_ab[1],xx_opti)

plt.figure()
plt.hold(True)
plt.plot(tt[Tstart:Tstart+4*OSR],uu_filt[Tstart:Tstart+4*OSR],'b',\
    label="Filtered input")
plt.plot(tt[Tstart:Tstart+4*OSR],xx_opti_filt[Tstart:Tstart+4*OSR],'r',\
    label="Filtered Optimal DSM output")
plt.legend(loc="best")
plt.suptitle("Portion of time domain behavior")

print("Filtering input")
uud_filt=sp.signal.lfilter(hz_ab[0],hz_ab[1],uud)

noise_power_opti_2=np.sum((uud_filt[Tstart:]-xx_opti_filt[Tstart:])**2)/ \
    (Tstop-Tstart)

print("Observed optimal noise level {} ({} dB)\nObserved SNR {} dB)".format( \
    noise_power_opti_2, dbp(noise_power_opti_2), \
        dbp(0.5*A1**2+0.5*A2**2)-dbp(noise_power_opti_2)))

plt.show()