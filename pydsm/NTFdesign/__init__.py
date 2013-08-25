# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
A collection of routines for the synthesis of Delta Sigma modulators.
=====================================================================

The delsig_optimize_NTF variable specifies whether to use optimization
techniques in the delsig_synthesizeNTF function.
"""

from . import filter_based
from . import weighting
from . import delsig
from . import psychoacoustic