# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

"""
A collection of routines for the synthesis of Delta Sigma modulators.
=====================================================================

The delsig_optimize_NTF variable specifies whether to use optimization 
techniques in the delsig_synthesizeNTF function.
"""

# Promote some functions/global variables to the synthesis namespace
from ..delsig import synthesizeNTF as delsig_synthesizeNTF
from ..delsig import optimize_NTF as delsig_optimize_NTF
import filter_based