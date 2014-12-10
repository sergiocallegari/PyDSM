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

u"""
Helper functions for NTF design (:mod:`pydsm.NTFdesign.helpers`)
================================================================

This modules provides helper functions for the design of the NTF
of ΔΣ modulators.

.. currentmodule:: pydsm.NTFdesign.helpers

.. autosummary::
   :toctree: generated/

   maxflat_fir_zeros     -- Zeros of a maxflat FIR transfer function
   spread_fir_uc_zeros   -- Zeros spread on unit circle according to cost
"""


from __future__ import division

import numpy as np
from scipy.optimize import minimize
from ..utilities import split_options, strip_options

__all__ = ["maxflat_fir_zeros", "spread_fir_uc_zeros"]


def maxflat_fir_zeros(order, alpha):
    """
    Compute the zeros of a maxflat FIR transfer function.

    The computed FIR transfer function is used in the DELSIG
    :func:`pydsm.delsig.synthesizeDSM()` design method
    (also known as :func:`pydsm.NTFdesign.ntf_schreier()`)
    for the denominator of the noise transfer function.

    Parameters
    ----------
    order : int
        the transfer function order
    alpha : float
        parameter controlling the steepness of the magnitude
        response

    Returns
    -------
    zeros : ndarray
        an array of complex roots of the FIR transfer function
    """
    x = 1./np.sqrt(alpha)
    me2 = -0.5*(x**(2./order))
    w = (2*np.arange(1, order+1)+1)*np.pi/order
    mb2 = 1+me2*np.exp(1j*w)
    p = mb2 - np.sqrt(mb2**2-1)
    # Reflect poles that fall out of the unit circle
    out = (np.abs(p) > 1)
    p[out] = 1/p[out]
    return p


def spread_fir_uc_zeros(order, OSR, cf, cf_args=[], cf_kwargs={}, **options):
    """
    Compute the best spreading of zerors on the unit circle.

    The computed FIR transfer function is optimal according to a
    criterion expressed by the cost function ``cf``.

     The computed FIR transfer function is used in the DELSIG
    :func:`pydsm.delsig.synthesizeDSM()` design method
    (also known as :func:`pydsm.NTFdesign.ntf_schreier()`)
    for the numerator of the noise transfer function.

    Parameters
    ----------
    order : int
        the transfer function order
    OSR : float
        the oversampling ratio
    cf : function
        cost function for the optimization. Takes a transfer
        function in zpk form as the first argument plus more
        arguments as required
    cf_args: list
        positional args of function ``cf``
    cf_kwargs: dict
        keyword args of function ``cf``

    Returns
    -------
    zeros : ndarray
        an array of complex roots of the FIR transfer function

    Other parameters
    ----------------
    L-BFGS-B_xxx : various types, optional
        Parameters prefixed by ``L-BFGS-B_`` are passed to the ``F-BFGS-B``
        optimizer. Allowed options are:

        ``L_BFGS_B_ftol``
            stop condition for the minimization
        ``L_BFGS_B_gtol``
            gradient stop condition for the minimization
        ``L_BFGS_B_maxcor``
            max number of variables used in hessian approximation
        ``L_BFGS_B_maxiter``
            max number of iterations
        ``L_BFGS_B_maxfun``
            max number of function evaluations
        ``L_BFGS_B_eps``
            Step size used for numerical approximation of the jacobian

        Do not use other options since they could break the minimizer in
        unexpected ways. Defaults can be set by changing the function
        ``default_options`` attribute.

    Notes
    -----
    The system is implicitly assumed to be low-pass. Hence, the zeros
    are spread on the unit circle in the [0, pi/OSR] range.

    See Also
    --------
    scipy.optimize.minimize :  Internally used minimizer
    """
    def dof2zeros(xx):
        zeros[0:xl] = np.exp(1j*xx)
        zeros[xl:2*xl] = zeros[0:xl].conj()
        if order % 2 is 1:
            zeros[-1] = 1
        return zeros

    def mf(xx):
        return np.log10(cf((dof2zeros(xx), np.zeros(order), 1),
                           *cf_args, **cf_kwargs))

    # Manage optional parameters
    opts = spread_fir_uc_zeros.default_options.copy()
    opts.update(options)
    o = split_options(opts, ['L_BFGS_B_'])
    lbfgsb_opts = strip_options(o, 'L_BFGS_B_')

    xl = order // 2
    zeros = np.zeros(order, dtype=complex)
    xx = minimize(mf, np.linspace(np.pi/OSR/order, np.pi/OSR, xl),
                  method='l-bfgs-b', options=lbfgsb_opts,
                  bounds=[(0, np.pi/OSR)] * xl).x
    return dof2zeros(xx)

spread_fir_uc_zeros.default_options = {"L_BFGS_B_ftol": 2.220446049250313e-09,
                                       "L_BFGS_B_gtol": 1e-05,
                                       "L_BFGS_B_maxcor": 10,
                                       "L_BFGS_B_maxiter": 15000,
                                       "L_BFGS_B_maxfun": 15000,
                                       "L_BFGS_B_eps": 1E-8}
