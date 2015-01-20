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

from __future__ import division, print_function

from numpy.testing import TestCase, run_module_suite
import numpy as np
from pydsm.NTFdesign import ntf_hybrid_weighting
from pydsm.relab import cplxpair
from nose.plugins.skip import SkipTest

__all__ = ["TestNTF_Hybrid"]


class TestNTF_Hybrid(TestCase):

    def setUp(self):
        # This test emulates a Schreier-type design using the hybrid
        # design method

        # Set the main design parameters
        self.order = 3
        self.OSR = 64

        # Set the NTF z, p, k that would be returned by Scheier's method
        self.e_k = 1
        e_z = [1.0000, 0.9993 - 0.0382j, 0.9993 + 0.0382j]
        self.e_z = cplxpair(e_z)
        e_p = [0.6692, 0.7652 - 0.2795j, 0.7652 + 0.2795j]
        self.e_p = cplxpair(e_p)

    # Prepare the weighting function for the hybrid method
    def w(self, f):
        return 1. if f <= 0.5/self.OSR else 1E-12

    def test_ntf_hybrid_tinoco(self):
        try:
            import cvxpy_tinoco     # analysis:ignore
        except:
            raise SkipTest("Modeler 'cvxpy_old' not installed")
        z, p, k = ntf_hybrid_weighting(self.order, self.w, H_inf=1.5,
                                       poles=self.e_p,
                                       show_progress=False,
                                       modeler='cvxpy_old',
                                       cvxopt_opts={"reltol": 1E-14,
                                                    "abstol": 1E-16})
        np.testing.assert_allclose(k, self.e_k, 1e-6)
        np.testing.assert_allclose(z, self.e_z, 3e-4)
        np.testing.assert_allclose(p, self.e_p, 3e-4)

    def test_ntf_hybrid_cvxpy(self):
        try:
            import cvxpy     # analysis:ignore
        except:
            raise SkipTest("Modeler 'cvxpy' not installed")
        z, p, k = ntf_hybrid_weighting(self.order, self.w, H_inf=1.5,
                                       poles=self.e_p,
                                       show_progress=False,
                                       modeler='cvxpy',
                                       cvxopt_opts={"reltol": 1E-14,
                                                    "abstol": 1E-16})
        np.testing.assert_allclose(k, self.e_k, 1e-6)
        np.testing.assert_allclose(z, self.e_z, 3e-4)
        np.testing.assert_allclose(p, self.e_p, 3e-4)

    def test_ntf_hybrid_picos(self):
        try:
            import picos     # analysis:ignore
        except:
            raise SkipTest("Modeler 'picos' not installed")
        z, p, k = ntf_hybrid_weighting(self.order, self.w, H_inf=1.5,
                                       poles=self.e_p,
                                       show_progress=False,
                                       modeler='picos',
                                       cvxopt_opts={"reltol": 1E-14,
                                                    "abstol": 1E-16})
        np.testing.assert_allclose(k, self.e_k, 1e-6)
        np.testing.assert_allclose(z, self.e_z, 3e-4)
        np.testing.assert_allclose(p, self.e_p, 3e-4)


if __name__ == '__main__':
    run_module_suite()
