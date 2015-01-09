# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
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

"""
Exceptions and warnings (:mod:`pydsm.exceptions`)
=================================================

This module defines the custom exceptions and warnings used in PyDSM

.. currentmodule:: pydsm.exceptions

Warning classes
---------------

.. autosummary::
   :toctree: generated/

    PyDsmWarning  -- base warning class
    PyDsmApproximationWarning  -- Approximation warning
    PyDsmSlowPathWarning  -- Slow path warning
    PyDsmPendingDeprecationWarning  -- Pending deprecation warning
    PyDsmDeprecationWarning  -- Deprecation warning


Error classes
-------------

No error classes are currently defined. Standard Python errors are risen.
"""


class PyDsmWarning(UserWarning):
    """Base class of warnings in PyDsm
    """
    pass


class PyDsmApproximationWarning(PyDsmWarning):
    """Class to warn about approximation in computation
    """
    pass


class PyDsmSlowPathWarning(PyDsmWarning):
    """Class to warn about slow code paths
    """
    pass


class PyDsmPendingDeprecationWarning(PyDsmWarning, PendingDeprecationWarning):
    """Class to warn about features pending deprecation
    """
    pass


class PyDsmDeprecationWarning(PyDsmWarning, DeprecationWarning):
    """Class to warn about deprecated features
    """
    pass
