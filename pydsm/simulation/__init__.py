# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
A collection of routines for the simulation of Delta Sigma modulators.
======================================================================

The delsig_use_fast_simulator variable specifies whether to use faster cython
code for the delsig_simulateDSM modulator simulator.
"""

# Promote some functions/global variables to the simulation namespace
from ..delsig import simulateDSM
from ..delsig import ds_quantize
from ..delsig import use_fast_simulator
