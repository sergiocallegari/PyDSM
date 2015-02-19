# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari, Federico Bizzarri
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

"""Demo code for a delta sigma modulator driving an electric motor at
different mechanical loads (and slips).

Copyright (c) Sergio Callegari, Federico Bizzarri 2012
All rights reserved.
"""

# Note. This code returns numeric results that are not fully consistent with
# those in the ICECS 2012 paper. Specifically, the SNR values
# evaluated via the NTF on the paper are better than those obtained by
# this code by a factor 2. In previous versions of the code, the results
# where identical due to a missing factor 2 in the quantization noise gain
# evaluation. Probably, this error slipped in because it makes the NTF-based
# results and the 'real' time-domain-simulation results more similar. In fact,
# the 'real' modulator behavior is better than predicted by the NTF model
# because in the specific example the quantization noise is not exactly white,
# but slightly blue (for some reason it gets more power at higher frequencies).


from __future__ import division, print_function

import numpy as np
import scipy as sp
__import__("scipy.signal")
import matplotlib.pyplot as plt
from pydsm.delsig import dbv, dbp, evalTF, synthesizeNTF, simulateDSM
from pydsm.ir import impulse_response
from pydsm.NTFdesign import quantization_noise_gain
from pydsm.NTFdesign.legacy import q0_from_filter_ir
from pydsm.NTFdesign.weighting import ntf_fir_from_q0


# Linearized motor model and default parameters
class Motor:
    P = 4       # Number of motor poles
    B = 25E-3   # Damping constant (dissipation due to windage and friction)
    J = 25E-3   # Moment of inertia
    Rs = 17.7   # Stator dissipative effects
    Rr = 13.8   # Rotor dissipative effects
    Ls = 459.2E-3   # Total 3-phase stator inductance
    Lr = 457.0E-3   # Total 3-phase rotor inductance
    Lm = 442.5E-3   # Magnetizing inductance
    fw = 50     # Nominal driving voltage angular frequency
    Vw = 320    # Nominal driving voltage peak amplitude

    def g0(self, sigma):
        return 1/self.Rs

    def tf_s(self, sigma, normalize=False):
        # Returns the linearized dynamical model of the motor
        # as a filter
        # sigma is slip coefficient
        b = np.array([sigma/self.Ls,
                      self.Rr/(self.Ls*self.Lr)])
        a = np.array([sigma*(1-self.Lm**2/(self.Ls*self.Lr)),
                      self.Rr/self.Lr+sigma*self.Rs/self.Ls,
                      self.Rr*self.Rs/(self.Lr*self.Ls)])
        if normalize:
            return (b/self.g0(sigma), a)
        else:
            return (b, a)

    def tf_z(self, sigma, fs, normalize=False):
        # Returns the linearized dynamical model of the motor
        # as a discrete time filter
        # sigma is slip coefficient
        # fs is the sample rate
        (btc, atc) = self.tf_s(sigma, normalize)
        return sp.signal.bilinear(btc, atc, fs)

motor = Motor()

# Motor data and Delta Sigma setup
fmot = 50.         # Nominal frequency of electrical drive
OSR = 1000         # Oversampling ratio
fphi = fmot*2*OSR  # Clock of Digital Delta Sigma Modulator
A = 190./320       # Normalized amplitude of the electric drive

# Approximation level (order) for the NTF
P = 8
# ... the same for the Delsig based NTF
DELSIG_P = 4

# Parameters to check the operation of the modulator in the time domain
Tstop = 2*OSR*100
Tstart = 2*OSR*80
# Intervals to make example plots of time domain data
T1p = Tstart
T2p = Tstart+4*OSR
dither_sigma = 1e-3

# The motor transfer functions to check (at different load/slip)
hs0043 = motor.tf_s(0.043)
hs02 = motor.tf_s(0.2)
hs06 = motor.tf_s(0.6)
hz0043 = motor.tf_z(0.043, fphi, True)
hz02 = motor.tf_z(0.2, fphi, True)
hz06 = motor.tf_z(0.6, fphi, True)

# Frequency axis for frequency domain plots (log scale)
fmin_log10 = np.ceil(np.log10(fmot/OSR))
fmax_log10 = np.log10(fphi/2)
ff = np.logspace(fmin_log10, fmax_log10, (fmax_log10-fmin_log10)*100)

# Plot frequency response of the motor
yyz0043 = evalTF(hz0043, np.exp(2j*np.pi*ff/fphi))
yyz02 = evalTF(hz02, np.exp(2j*np.pi*ff/fphi))
yyz06 = evalTF(hz06, np.exp(2j*np.pi*ff/fphi))

plt.figure()
plt.semilogx(ff, dbv(np.abs(yyz0043)), 'b', label='$\sigma=0.043$')
plt.semilogx(ff, dbv(np.abs(yyz02)), 'r', label='$\sigma=0.2$')
plt.semilogx(ff, dbv(np.abs(yyz06)), 'g', label='$\sigma=0.6$')
plt.xlim(xmax=10**fmax_log10)
plt.ylim(ymin=-50)
plt.xlabel('$f$ [Hz]')
plt.ylabel('[dB]')
plt.suptitle("Magnitude response of motor linearized model")
plt.legend()
plt.ion()
plt.show()
plt.ioff()

# Compute optimal and benchmark NTFs
hz0043_ir = impulse_response(hz0043)
hz02_ir = impulse_response(hz02)
hz06_ir = impulse_response(hz06)
q0_0043 = q0_from_filter_ir(P, hz0043_ir)
ntf0043 = ntf_fir_from_q0(q0_0043)
q0_02 = q0_from_filter_ir(P, hz02_ir)
ntf02 = ntf_fir_from_q0(q0_02)
q0_06 = q0_from_filter_ir(P, hz06_ir)
ntf06 = ntf_fir_from_q0(q0_06)
delsig_ntf = synthesizeNTF(DELSIG_P, OSR, 3, 1.5, 0)

# Plot the NTFs
nyz0043 = evalTF(ntf0043, np.exp(2j*np.pi*ff/fphi))
nyz02 = evalTF(ntf02, np.exp(2j*np.pi*ff/fphi))
nyz06 = evalTF(ntf06, np.exp(2j*np.pi*ff/fphi))
delsig_ny = evalTF(delsig_ntf, np.exp(2j*np.pi*ff/fphi))

plt.figure()
plt.semilogx(ff, dbv(np.abs(nyz0043)), 'b',
             label='opt @ $\sigma=0.043$')
plt.semilogx(ff, dbv(np.abs(nyz02)), 'r',
             label='opt @ $\sigma=0.2$')
plt.semilogx(ff, dbv(np.abs(nyz06)), 'g',
             label='opt @ $\sigma=0.6$')
plt.semilogx(ff, dbv(np.abs(delsig_ny)), 'y', label='benchmark')
plt.xlabel('$f$ [Hz]')
plt.ylabel('[dB]')
plt.suptitle("NTF magnitude response")
plt.xlim(xmax=10**fmax_log10)
plt.ylim(-225, 5)
plt.legend(loc=4)

print("Computing the SNR in the various cases")
noise_power = 1./3  # Quantization noise is power is Delta^2/12
# low slip
n_pow_0043_0043 = noise_power*quantization_noise_gain(hz0043, ntf0043)
n_pow_0043_02 = noise_power*quantization_noise_gain(hz0043, ntf02)
n_pow_0043_06 = noise_power*quantization_noise_gain(hz0043, ntf06)
n_pow_0043_delsig = noise_power*quantization_noise_gain(hz0043, delsig_ntf)
# mid slip
n_pow_02_0043 = noise_power*quantization_noise_gain(hz02, ntf0043)
n_pow_02_02 = noise_power*quantization_noise_gain(hz02, ntf02)
n_pow_02_06 = noise_power*quantization_noise_gain(hz02, ntf06)
n_pow_02_delsig = noise_power*quantization_noise_gain(hz02, delsig_ntf)
# high slip
n_pow_06_0043 = noise_power*quantization_noise_gain(hz06, ntf0043)
n_pow_06_02 = noise_power*quantization_noise_gain(hz06, ntf02)
n_pow_06_06 = noise_power*quantization_noise_gain(hz06, ntf06)
n_pow_06_delsig = noise_power*quantization_noise_gain(hz06, delsig_ntf)

sig_pow = A**2/2  # Input signal power
# Output signal power at 3 slip conditions
sig_pow_0043 = sig_pow*np.abs(evalTF(hz0043, np.exp(2j*np.pi*fmot/fphi)))**2
sig_pow_02 = sig_pow*np.abs(evalTF(hz02, np.exp(2j*np.pi*fmot/fphi)))**2
sig_pow_06 = sig_pow*np.abs(evalTF(hz06, np.exp(2j*np.pi*fmot/fphi)))**2

print("Noise power and SNR for motor at slip=0.043")
print("ntf0043\tntf02\tntf06\tdelsig")
print("{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}".format(
    n_pow_0043_0043, n_pow_0043_02, n_pow_0043_06, n_pow_0043_delsig))
print("{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB".format(
    dbp(sig_pow_0043/n_pow_0043_0043), dbp(sig_pow_0043/n_pow_0043_02),
    dbp(sig_pow_0043/n_pow_0043_06), dbp(sig_pow_0043/n_pow_0043_delsig)))
print("Noise power and SNR for motor at slip=0.2")
print("ntf0043\tntf02\tntf06\tdelsig")
print("{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}".format(
    n_pow_02_0043, n_pow_02_02, n_pow_02_06, n_pow_02_delsig))
print("{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB".format(
    dbp(sig_pow_02/n_pow_02_0043), dbp(sig_pow_02/n_pow_02_02),
    dbp(sig_pow_02/n_pow_02_06), dbp(sig_pow_02/n_pow_02_delsig)))
print("Noise power and SNR for motor at slip=0.6")
print("ntf0043\tntf02\tntf06\tdelsig")
print("{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}".format(
    n_pow_06_0043, n_pow_06_02, n_pow_06_06, n_pow_06_delsig))
print("{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB".format(
    dbp(sig_pow_06/n_pow_06_0043), dbp(sig_pow_06/n_pow_06_02),
    dbp(sig_pow_06/n_pow_06_06), dbp(sig_pow_06/n_pow_06_delsig)))

# Repeating the tests in the time domain with the nonlinear modulator model
# and the linearized motor model
# Setting up DS simulation
tt = np.asarray(range(int(Tstop)))
uu = A*np.sin(2*np.pi*fmot/fphi*tt)
dither = np.random.randn(len(uu))*dither_sigma
uud = uu+dither
# Run the simulations
print("Simulating modulator with NTF optimized for sigma=0.043")
xx_opti0043 = simulateDSM(uud, ntf0043)[0]
print("Simulating modulator with NTF optimized for sigma=0.2")
xx_opti02 = simulateDSM(uud, ntf02)[0]
print("Simulating modulator with NTF optimized for sigma=0.6")
xx_opti06 = simulateDSM(uud, ntf06)[0]
print("Simulating modulator with benchmark NTF")
xx_delsig = simulateDSM(uud, delsig_ntf)[0]
print("Applying the filters from the linearized motor model")
print("Using the transfer function corresponding to sigma=0.043")
uud_filt0043 = sp.signal.lfilter(hz0043[0], hz0043[1], uud)
xx_opti_filt0043_0043 = sp.signal.lfilter(hz0043[0], hz0043[1], xx_opti0043)
xx_opti_filt0043_02 = sp.signal.lfilter(hz0043[0], hz0043[1], xx_opti02)
xx_opti_filt0043_06 = sp.signal.lfilter(hz0043[0], hz0043[1], xx_opti06)
xx_delsig_filt0043 = sp.signal.lfilter(hz0043[0], hz0043[1], xx_delsig)
print("Using the filter for sigma=0.2")
uud_filt02 = sp.signal.lfilter(hz02[0], hz02[1], uud)
xx_opti_filt02_0043 = sp.signal.lfilter(hz02[0], hz02[1], xx_opti0043)
xx_opti_filt02_02 = sp.signal.lfilter(hz02[0], hz02[1], xx_opti02)
xx_opti_filt02_06 = sp.signal.lfilter(hz02[0], hz02[1], xx_opti06)
xx_delsig_filt02 = sp.signal.lfilter(hz02[0], hz02[1], xx_delsig)
print("Using the filter for sigma=0.6")
uud_filt06 = sp.signal.lfilter(hz06[0], hz06[1], uud)
xx_opti_filt06_0043 = sp.signal.lfilter(hz06[0], hz06[1], xx_opti0043)
xx_opti_filt06_02 = sp.signal.lfilter(hz06[0], hz06[1], xx_opti02)
xx_opti_filt06_06 = sp.signal.lfilter(hz06[0], hz06[1], xx_opti06)
xx_delsig_filt06 = sp.signal.lfilter(hz06[0], hz06[1], xx_delsig)
print("Computing the merit factors")
noise_power_opti_0043_0043 = np.sum(
    (uud_filt0043[Tstart:]-xx_opti_filt0043_0043[Tstart:])**2)/(Tstop-Tstart)
noise_power_opti_0043_02 = np.sum(
    (uud_filt0043[Tstart:]-xx_opti_filt0043_02[Tstart:])**2)/(Tstop-Tstart)
noise_power_opti_0043_06 = np.sum(
    (uud_filt0043[Tstart:]-xx_opti_filt0043_06[Tstart:])**2)/(Tstop-Tstart)
noise_power_delsig_0043 = np.sum(
    (uud_filt0043[Tstart:]-xx_delsig_filt0043[Tstart:])**2)/(Tstop-Tstart)
noise_power_opti_02_0043 = np.sum(
    (uud_filt02[Tstart:]-xx_opti_filt02_0043[Tstart:])**2)/(Tstop-Tstart)
noise_power_opti_02_02 = np.sum(
    (uud_filt02[Tstart:]-xx_opti_filt02_02[Tstart:])**2)/(Tstop-Tstart)
noise_power_opti_02_06 = np.sum(
    (uud_filt02[Tstart:]-xx_opti_filt02_06[Tstart:])**2)/(Tstop-Tstart)
noise_power_delsig_02 = np.sum(
    (uud_filt02[Tstart:]-xx_delsig_filt02[Tstart:])**2)/(Tstop-Tstart)
noise_power_opti_06_0043 = np.sum(
    (uud_filt06[Tstart:]-xx_opti_filt06_0043[Tstart:])**2)/(Tstop-Tstart)
noise_power_opti_06_02 = np.sum(
    (uud_filt06[Tstart:]-xx_opti_filt06_02[Tstart:])**2)/(Tstop-Tstart)
noise_power_opti_06_06 = np.sum(
    (uud_filt06[Tstart:]-xx_opti_filt06_06[Tstart:])**2)/(Tstop-Tstart)
noise_power_delsig_06 = np.sum(
    (uud_filt06[Tstart:]-xx_delsig_filt06[Tstart:])**2)/(Tstop-Tstart)
sig_powx_0043 = np.sum(uud_filt0043[Tstart:]**2)/(Tstop-Tstart)
sig_powx_02 = np.sum(uud_filt02[Tstart:]**2)/(Tstop-Tstart)
sig_powx_06 = np.sum(uud_filt06[Tstart:]**2)/(Tstop-Tstart)
print("Noise power and SNR for motor at slip=0.043")
print("ntf0043\tntf02\tntf06\tdelsig")
print("{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}".format(
    noise_power_opti_0043_0043, noise_power_opti_0043_02,
    noise_power_opti_0043_06, noise_power_delsig_0043))
print("{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB".format(
    dbp(sig_powx_0043/noise_power_opti_0043_0043),
    dbp(sig_powx_0043/noise_power_opti_0043_02),
    dbp(sig_powx_0043/noise_power_opti_0043_06),
    dbp(sig_powx_0043/noise_power_delsig_0043)))
print("Noise power and SNR for motor at slip=0.2")
print("ntf0043\tntf02\tntf06\tdelsig")
print("{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}".format(
    noise_power_opti_02_0043, noise_power_opti_02_02,
    noise_power_opti_02_06, noise_power_delsig_02))
print("{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB".format(
    dbp(sig_powx_02/noise_power_opti_02_0043),
    dbp(sig_powx_02/noise_power_opti_02_02),
    dbp(sig_powx_02/noise_power_opti_02_06),
    dbp(sig_powx_02/noise_power_delsig_02)))
print("Noise power and SNR for motor at slip=0.6")
print("ntf0043\tntf02\tntf06\tdelsig")
print("{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}".format(
    noise_power_opti_06_0043, noise_power_opti_06_02,
    noise_power_opti_06_06, noise_power_delsig_06))
print("{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB\t{:2.2f}dB".format(
    dbp(sig_powx_06/noise_power_opti_06_0043),
    dbp(sig_powx_06/noise_power_opti_06_02),
    dbp(sig_powx_06/noise_power_opti_06_06),
    dbp(sig_powx_06/noise_power_delsig_06)))

# Plot portions of the time domain behavior
plt.figure()
plt.suptitle("Current profile at sigma=0043")
plt.plot(tt[T1p:T2p], xx_opti_filt0043_0043[T1p:T2p], 'b', label='optimal NTF')
plt.plot(tt[T1p:T2p], xx_delsig_filt0043[T1p:T2p], 'r', label='benchmark NTF')
plt.plot(tt[T1p:T2p], uud_filt0043[T1p:T2p], 'g', label='ideal behavior')
plt.legend()
plt.xlabel('time')
plt.ylabel('normalized currents')

plt.show()
