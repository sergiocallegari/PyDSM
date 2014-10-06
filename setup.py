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
import os
from setuptools import setup, Command, Extension, find_packages
from Cython.Distutils import build_ext
import platform
import numpy as np

__version__ = ''
exec(open('pydsm/_version.py').read())


def read_from_here(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as fp:
        return fp.read()


class test (Command):
    description = "Test the pydsm distribution prior to install"

    user_options = [
        ('test-file=', None,
         'Testfile to run in the test directory'),
        ]

    def initialize_options(self):
        self.build_base = 'build'
        self.test_dir = 'test'
        self.test_file = 'test_all'

    def finalize_options(self):
        build = self.get_finalized_command('build')
        self.build_purelib = build.build_purelib
        self.build_platlib = build.build_platlib

    def run(self):
        # Invoke the 'build' command
        self.run_command('build')
        # remember old sys.path to restore it afterwards
        old_path = sys.path[:]
        # extend sys.path
        sys.path.insert(0, self.build_purelib)
        sys.path.insert(0, self.build_platlib)
        sys.path.insert(0, self.test_dir)
        # build include path for test
        TEST = __import__(self.test_file)
        suite = TEST.unittest.TestLoader().loadTestsFromModule(TEST)
        TEST.unittest.TextTestRunner(verbosity=2).run(suite)
        sys.path = old_path[:]

# Prepare the extension modules
ext_modules = [
    Extension('pydsm.delsig._simulateDSM_cblas',
              ['pydsm/delsig/_simulateDSM_cblas.pyx'],
              include_dirs=[np.get_include()],
              libraries=['cblas']),
    Extension('pydsm.delsig._simulateDSM_scipy_blas',
              ['pydsm/delsig/_simulateDSM_scipy_blas.pyx'],
              include_dirs=[np.get_include()])]

description = 'Python Based Delta-Sigma modulator design tools'
# Long description can contain restructured text and goes on display
# on Pypi
long_description = (read_from_here('README') +
                    '\n\n' +
                    read_from_here('CHANGELOG'))

# Special requirements for the windows platform
if platform.system() != 'Linux':
    # The cblas simulator is only built in Linux
    ext_modules = [
        Extension('pydsm.delsig._simulateDSM_scipy_blas',
                  ['pydsm/delsig/_simulateDSM_scipy_blas.pyx'],
                  include_dirs=[np.get_include()])]

setup(
    name='pydsm',
    version=__version__,
    description=description,
    long_description=long_description,
    author='Sergio Callegari',
    author_email='sergio.callegari@unibo.it',
    url='http://pydsm.googlecode.com',
    license='GNU General Public License v3 or later (GPLv3+)',
    platforms=['Linux', 'Windows', 'Mac'],
    packages=find_packages(),
    ext_modules=ext_modules,
    requires=['scipy(>=0.10.1)',
              'numpy(>=1.6.1)',
              'matplotlib(>= 1.1.0)',
              'cvxopt(>=1.1.4)',
              'cython(>=0.16)'],
    cmdclass={'test': test,
              'build_ext': build_ext},
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        ('License :: OSI Approved :: '
         'GNU General Public License v3 or later (GPLv3+)'),
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Cython',
        'Topic :: Education',
        ('Topic :: Scientific/Engineering :: '
         'Electronic Design Automation (EDA)')
        ]
)
