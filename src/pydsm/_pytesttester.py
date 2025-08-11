# Copyright © 2024, Sergio Callegari
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


# Based on code Copyright (c) 2005-2024, NumPy Developers.
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:

#    * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#    * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#    * Neither the name of the NumPy Developers nor the names of any
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
Pytest test runner.

This module implements the ``test()`` function for PyDSM modules and is
based on the NumPy module with the same name

Warnings filtering and other runtime settings should be dealt with in the
``pytest.ini`` file in the repo root. The behavior of the test depends on
whether or not that file is found as follows:

* ``pytest.ini`` is present (develop mode)
    All warnings except those explicitly filtered out are raised as error.
* ``pytest.ini`` is absent (release mode)
    DeprecationWarnings and PendingDeprecationWarnings are ignored, other
    warnings are passed through.

In practice, tests run from the repo are run in development mode from an
inplace build with ``pytest pydsm``.

This module is imported by every subpackage.
"""
import sys
import os

__all__ = ['PytestTester']


def _show_pydsm_info():
    import pydsm
    print(f"PyDSM version {pydsm.__version__}")

class PytestTester:
    """
    Pytest test runner.

    Calling this test function finds and runs all tests associated with the
    module and all its sub-modules.

    Attributes
    ----------
    module_name : str
        The name of the (sub)package/module to test.

    Parameters
    ----------
    module_name : module name
        The name of the (sub)package/module to test.

    Notes
    -----
    This class is not publicly exposed.

    """

    def __init__(self, module_name, mode='tests'):
        self.module_name = module_name
        assert mode in ['tests', 'benchmarks']
        self.mode = mode

    def __call__(self, label='fast', verbose=1, extra_argv=None,
                 doctests=False, coverage=False, durations=-1, tests=None):
        """
        Run tests for module using pytest.

        Parameters
        ----------
        label : {'fast', 'full'}, optional
            Identifies the tests to run. When set to 'fast', tests decorated
            with `pytest.mark.slow` are skipped, when 'full', the slow marker
            is ignored.
        verbose : int, optional
            Verbosity value for test outputs, in the range 1-3. Default is 1.
        extra_argv : list, optional
            List with any extra arguments to pass to pytests.
        doctests : bool, optional
            Perform doctests
        coverage : bool, optional
            If True, report coverage. Default is False.
            Requires installation of (pip) pytest-cov.
        durations : int, optional
            If < 0, do nothing, If 0, report time of all tests, if > 0,
            report the time of the slowest `timer` tests. Default is -1.
        tests : test or list of tests
            Tests to be executed with pytest '--pyargs'

        Returns
        -------
        result : bool
            Return True on success, false otherwise.

        Notes
        -----
        Each PyDSM module exposes `test` in its namespace to run all tests for
        it. E.g. `pydsm.NTFdesign.test()`.
        """

        import pytest
        import warnings

        module = sys.modules[self.module_name]
        module_path = os.path.abspath(module.__path__[0])

        # setup the pytest arguments
        pytest_args = ["-l"]

        # offset verbosity. The "-q" cancels a "-v".
        pytest_args += ["-q"]

        if self.mode == 'tests':
            pytest_args += ["--ignore-glob", "**/benchmark_*.py"]
        else:
            pytest_args += ["--ignore-glob", "**/test_*.py"]

        if doctests:
            pytest_args += ["--doctest-modules"]

        if extra_argv:
            pytest_args += list(extra_argv)

        if verbose > 1:
            pytest_args += ["-" + "v"*(verbose - 1)]

        if coverage:
            pytest_args += ["--cov=" + module_path]

        if label == "fast":
            pytest_args += ["-m", "not slow"]

        elif label != "full":
            pytest_args += ["-m", label]

        if durations >= 0:
            pytest_args += [f"--durations={durations}"]

        if tests is None:
            tests = [self.module_name]

        pytest_args += ["--pyargs"] + list(tests)

        # run tests.
        _show_pydsm_info()

        try:
            code = pytest.main(pytest_args)
        except SystemExit as exc:
            code = exc.code

        return code == 0
