#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014–2025, Sergio Callegari
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


import sys
from setuptools import setup, Extension
from Cython.Build import cythonize
import platform
import numpy as np

if sys.version_info[:2] < (3, 10):
    raise RuntimeError("Python 3 supported for versions >= 3.10")

# Prepare the extension modules
ext_modules = cythonize([
    Extension(name='pydsm.delsig._simulateDSM_cblas',
              sources=['src/pydsm/delsig/_simulateDSM_cblas.pyx'],
              include_dirs=[np.get_include()]+["/usr/include/openblas/"],
              libraries=['blas'],
              define_macros=[('NPY_NO_DEPRECATED_API',
                              'NPY_1_7_API_VERSION')]),
    Extension(name='pydsm.delsig._simulateDSM_scipy_blas',
              sources=['src/pydsm/delsig/_simulateDSM_scipy_blas.pyx'],
              include_dirs=[np.get_include()],
              define_macros=[('NPY_NO_DEPRECATED_API',
                              'NPY_1_7_API_VERSION')])],
    compiler_directives={'language_level' : "3"})

# Special requirements for the windows platform
if platform.system() != 'Linux':
    # The cblas simulator is only built in Linux
    ext_modules = [
        Extension(name='pydsm.delsig._simulateDSM_scipy_blas',
                  sources=['src/pydsm/delsig/_simulateDSM_scipy_blas.pyx'],
                  include_dirs=[np.get_include()],
                  define_macros=[('NPY_NO_DEPRECATED_API',
                                  'NPY_1_7_API_VERSION')])]

setup(
    name='pydsm',
    setup_requires=["setuptools_scm"],
    ext_modules=ext_modules,
)
