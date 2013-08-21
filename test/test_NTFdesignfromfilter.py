# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

import unittest
import numpy as np
import scipy as sp
__import__("scipy.signal")
from pydsm.ir import impulse_response
from pydsm.delsig import evalTF
from pydsm.NTFdesign.filter_based import synthesize_ntf_from_filter_imp,\
    synthesize_ntf_from_filter_mag, quantization_noise_gain

__all__=["TestNTF_Filter"]

class TestNTF_Filter(unittest.TestCase):

    def setUp(self):
        pass

    def test_ntf_butt_bp8(self):
        # Generate filter.
        # 8th order bandpass filter
        # Freq. passed to butterworth is normalized between 0 and 1
        # where 1 is the Nyquist frequency
        fsig=1000.
        B=400.
        OSR=64
        fphi=B*OSR*2
        w0=2*fsig/fphi
        B0=2*B/fphi
        w1=(np.sqrt(B0**2+4*w0**2)-B0)/2
        w2=(np.sqrt(B0**2+4*w0**2)+B0)/2
        hz=sp.signal.butter(4, [w1,w2], 'bandpass', output='zpk')
        # Order
        order=12
        # Compute q0 in two ways
        ir=impulse_response(hz,db=80)
        def mr(f):
            return np.abs(evalTF(hz,np.exp(2j*np.pi*f)))
        ntf1=synthesize_ntf_from_filter_mag(order, mr,
                                            options={'show_progress':False})
        ntf2=synthesize_ntf_from_filter_imp(order, ir,
                                            options={'show_progress':False})
        mf1=quantization_noise_gain(ntf1,hz)
        mf2=quantization_noise_gain(ntf2,hz)
        np.testing.assert_almost_equal(mf1,mf2,decimal=12)

if __name__ == '__main__':
    unittest.main()
