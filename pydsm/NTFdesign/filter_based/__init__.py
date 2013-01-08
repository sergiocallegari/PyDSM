# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
Output filter based NTF synthesis
=================================

This modules provides code for the synthesis of the modulator NTF, based
on the filter placed after the modulator for the removal of the quantization
noise.

The algorithms used in these routines are extensively described in

Sergio Callegari, Federico Bizzarri "Output Filter Aware Optimization of the
Noise Shaping Properties of ΔΣ Modulators via Semi-Definite Programming",
IEEE Transactions on Circuits and Systems I: Regular Papers.
"""

from _q0_from_filter import *
from _quantization_noise_gain import *
from _ntf_from_filter import *