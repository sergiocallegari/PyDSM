# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
Library for the design and simulation of Delta-Sigma modulators
===============================================================

Copyright (c) 2012, Sergio Callegari
Some code ported from the DELSIG toolbox, copyrighted by R. Schreier.
"""

# Read version info
from .version import __version__

# Promote some key functions to the pydsm namespace
from . import delsig
from . import NTFdesign
from . import simulation
from . import audio_weightings
