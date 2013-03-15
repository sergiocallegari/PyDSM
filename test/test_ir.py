# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

# -*- coding: utf-8 -*-

import unittest
import numpy as np
import scipy as sp
__import__("scipy.signal")
from pydsm.ir import impulse_response
from pydsm.delsig import evalTF
from pydsm.utilities import db

__all__=["TestIR"]

class TestIR(unittest.TestCase):

    def setUp(self):
        pass

    def test_butt_bp8_ir(self):
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
        ir=impulse_response(hz,db=80)
        ff=np.logspace(np.log10(w0/4), np.log10(w0*4), 100)
        vv1=evalTF(hz,np.exp(2j*np.pi*ff))
        vv2=evalTF((ir,[1]),np.exp(2j*np.pi*ff))
        np.testing.assert_allclose(db(vv1/vv2),np.zeros_like(vv1),\
            rtol=0,atol=3)

    def test_fir_ir(self):
        fir=np.asarray([1.0,0.5,0.25,0.125])
        zpk=(np.roots(fir),np.zeros(4),1)
        ir=impulse_response(zpk,db=80)
        np.testing.assert_allclose(fir,ir)

if __name__ == '__main__':
    unittest.main()
