#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014, Sergio Callegari
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
from distutils.core import setup, Command
from distutils.extension import Extension
from Cython.Build import cythonize
import platform
import numpy as np


# Find version
__version__=''
execfile('pydsm/_version.py')


class test (Command):
    description = "Test the pydsm distribution prior to install"

    user_options = [
        ('test-file=', None,
         'Testfile to run in the test directory'),
        ]

    def initialize_options (self):
        self.build_base = 'build'
        self.test_dir = 'test'
        self.test_file = 'test_all'

    def finalize_options (self):
        build = self.get_finalized_command('build')
        self.build_purelib = build.build_purelib
        self.build_platlib = build.build_platlib

    def run (self):
        # Invoke the 'build' command
        self.run_command ('build')
        # remember old sys.path to restore it afterwards
        old_path = sys.path[:]
        # extend sys.path
        sys.path.insert(0, self.build_purelib)
        sys.path.insert(0, self.build_platlib)
        sys.path.insert(0, self.test_dir)
        # build include path for test
        TEST=__import__(self.test_file)
        suite = TEST.unittest.TestLoader().loadTestsFromModule(TEST)
        TEST.unittest.TextTestRunner(verbosity=2).run(suite)
        sys.path = old_path[:]

# Prepare the extension modules
ext_modules=[
    Extension('pydsm.delsig._simulateDSM_cblas',
              ['pydsm/delsig/_simulateDSM_cblas.pyx']),
    Extension('pydsm.delsig._simulateDSM_scipy_blas',
              ['pydsm/delsig/_simulateDSM_scipy_blas.pyx'])]

description='Python Based Delta-Sigma modulator design tools'
long_description="""Python Based Delta-Sigma modulator design tools.

Based on the algorithms in Callegari, Bizzarri 'Output Filter Aware
Optimization of the Noise Shaping Properties of Delta-Sigma Modulators via
Semi-Definite Programming', IEEE Transactions on Circuits and Systems I,
2013 and others.

Portion of code ported to python from the DELSIG toolbox by R. Schreier.
"""

# Special requirements for the windows platform
if platform.system()=='Windows':
    # In windows, the cblas simulator is not built
    ext_modules=[
        Extension('pydsm.delsig._simulateDSM_scipy_blas',
                  ['pydsm/delsig/_simulateDSM_scipy_blas.pyx'],
                  include_dirs=[np.get_include()])]

setup(name='pydsm',
      version=__version__,
      description=description,
      author='Sergio Callegari',
      author_email='sergio.callegari@unibo.it',
      url='http://pydsm.googlecode.com',
      packages = ['pydsm', 'pydsm.simulation', 'pydsm.NTFdesign',
                  'pydsm.NTFdesign.filter_based', 'pydsm.delsig',
                  'pydsm.NTFdesign.weighting', 'cvxpy_tinoco',
		  'cvxpy_tinoco.functions', 'cvxpy_tinoco.procedures',
		  'cvxpy_tinoco.sets'],
      ext_modules=cythonize(ext_modules),
      requires=['scipy(>=0.10.1)',
                'numpy(>=1.6.1)',
                'matplotlib(>= 1.1.0)',
                'cvxopt(>=1.1.4)',
                'cython(>=0.16)'],
      cmdclass = {'test': test},
      license = 'GNU GPL version 3 or any later version',
      platforms = ['Linux','Windows','Mac'],
      long_description = long_description)
