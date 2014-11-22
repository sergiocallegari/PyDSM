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

import unittest

from pydsm.delsig import synthesizeNTF
from pydsm.delsig import plotPZ

__all__ = ["TestPlotPZ"]


class TestPlotPZ(unittest.TestCase):

    def setUp(self):
        pass

    def test_default(self):
        ntf = synthesizeNTF()
        plotPZ(ntf)

if __name__ == '__main__':
    unittest.main()
