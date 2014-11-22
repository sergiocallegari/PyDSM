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

from numpy.testing import TestCase, run_module_suite
from matplotlib import pyplot as plt
from pydsm.delsig import synthesizeNTF
from pydsm.delsig import plotPZ

__all__ = ["TestPlotPZ"]


class TestPlotPZ(TestCase):

    def setUp(self):
        self.saved_backend = plt.get_backend()
        plt.switch_backend('Agg')

    def tearDown(self):
        plt.switch_backend(self.saved_backend)

    def test_default(self):
        ntf = synthesizeNTF()
        plotPZ(ntf)

if __name__ == '__main__':
    run_module_suite()
