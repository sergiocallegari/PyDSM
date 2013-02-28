#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

import sys
import os
import shutil
from distutils.core import setup
from distutils.core import Command
from distutils import ccompiler
from distutils.extension import Extension
from Cython.Build import cythonize
import platform
import numpy


# Find version
execfile('pydsm/version.py')


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


# We want the default compiler to be mingw32 in windows
ccompiler._default_compilers=\
        (('nt', 'mingw32'),)+ccompiler._default_compilers

# Prepare the extension modules
ext_modules=[
    Extension('pydsm.delsig._simulateDSM_fast',
              ['pydsm/delsig/_simulateDSM_fast.pyx'])]

# Fix stuff for Windows
plat=platform.system()
if plat=='Windows':
    ext_modules=[
        Extension('pydsm.delsig._simulateDSM_fast',
                  ['pydsm/delsig/_simulateDSM_fast.pyx'],
                  include_dirs=[numpy.get_include(), './ATLAS'],
                  library_dirs=['./ATLAS'],
                  libraries=['atlas'],
                  define_macros=[('__USE_MINGW_ANSI_STDIO','1')])]

setup(name='pydsm',
      version=__version__,
      description=u'Python Based ΔΣ modulator design tools',
      author='Sergio Callegari',
      author_email='sergio.callegari@unibo.it',
      url='http://code.google.com/p/pydsm',
      packages = ['pydsm', 'pydsm.simulation', 'pydsm.NTFdesign',
                  'pydsm.NTFdesign.filter_based', 'pydsm.delsig'],
      ext_modules=cythonize(ext_modules),
      requires=['scipy(>=0.10.1)',
                'numpy(>=1.6.1)',
                'matplotlib(>= 1.1.0)',
                'cvxopt(>=1.1.4)',
                'cython(>=0.16)'],
      cmdclass = {'test': test},
      license = 'Simplified BSD License',
      platforms = ['Linux','Windows','Mac'],
      long_description = u"""
      Python Based ΔΣ modulator design tools.

    Based on the algorithms in Callegari, Bizzarri 'Output Filter Aware
    Optimization of the Noise Shaping Properties of ΔΣ Modulators via
    Semi-Definite Programming', IEEE Transactions on Circuits and Systems
    I, 2013 and others.

    Portion of code ported to python from the DELSIG toolbox by R. Schreier.
    """)
