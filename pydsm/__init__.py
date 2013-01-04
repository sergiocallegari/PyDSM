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
from pydsm.version import __version__, __delsig_version__

# Promote some key functions to the pydsm namespace
import delsig
import NTFdesign
import simulation