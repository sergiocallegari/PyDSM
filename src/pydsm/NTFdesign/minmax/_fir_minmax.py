# -*- coding: utf-8 -*-

# Copyright (c) 2013, Sergio Callegari
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

import numpy as np
from warnings import warn
from ...exceptions import PyDsmDeprecationWarning
from ...utilities import digested_options

__all__ = ['ntf_fir_minmax', 'synthesize_ntf_minmax']


def ntf_fir_minmax(order=32, osr=32, H_inf=1.5, f0=0, zf=False,
                   **options):
    u"""Synthesize FIR NTF for LP or BP ΔΣ modulator by min-max optimization.

    The strategy aims at minimizing the peak value of the NTF in the signal
    band, while respecting the Lee criterion.

    Parameters
    ----------
    order : int, optional
        Order of the delta-sigma modulator. Defaults to 32.
    osr : real or array of reals, optional
        The oversampling ratio, based on the input signal bandwidth
        Defaults to 32. If there are multiple signal bands, multiple
        OSRs can be provided, one for each of them.
    H_inf : real, optional
        Max allowed peak value of the NTF. Used to enforce the Lee criterion.
        Defaults to 1.5.
    f0 : real, or array of reals, optional
        The normalized center frequency of the modulator. Value in [0,1/2].
        1/2 indicates half the sample frequency. Defaults to 0, indicating an
        LP modulator. An array of values can be provided for multiband design.
    zf : bool, optional
        Flag controlling the pre-assignement of NTF zeros. If ``False``, the
        design is practiced without any zero pre-assignment. If ``True``, a
        zero is pre-assigned at the modulator center-band. Defaults to False.

    Returns
    -------
    ntf : tuple
        noise transfer function in zpk form.

    cvxpy_opts : dictionary, optional
       A dictionary of options to use with the ``cvxpy`` modeling library.
       Allowed options include:

       ``override_kktsolver`` (bool)
           Whether to override the default ``cvxopt`` kkt solver using the
           ``chol`` kkt solver.
           Leave this at the default False setting, to avoid errors.
       ``solver`` (string)
           The solver backend to use. Either `cvxopt` or `scs`

    cvxopt_opts : dict, optional
        A dictionary of options for the ``cvxopt`` optimizer.
        Allowed options include:

        ``maxiters`` (int)
            Maximum number of iterations
        ``abstol`` (real)
            Absolute accuracy
        ``reltol`` (real)
            Relative accuracy
        ``feastol`` (real)
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break ``cvxopt`` in
        unexpected ways. These options can be passed when using the
        ``cvxpy_old`` modeler, the ``picos`` modeler or the ``cvxpy`` modeler
        with the ``cvxopt`` backend.
    scs_opts : dict, optional
        A dictionary of options for the ``scs`` optimizer.  Allowed options
        include:

        ``max_iters`` (int)
            Maximum number of iterations
        ``eps`` (real)
            Convergence tolerance
        ``alpha`` (real)
            Relaxation parameter
        ``normalize`` (bool)
            Whether to precondition data matrices
        ``use_indirect`` (bool)
            Whether to use indirect solver for KKT sytem (instead of direct)

       Do not use other options since they could break ``scs`` in
       unexpected ways. These options can be passed when using the
       ``cvxpy`` modeler with the ``scs`` backend.

    Notes
    -----
    The design strategy implemented in this module is described in the paper
    [1]_.

    Default values for the options not directly documented in the function
    call signature can be checked and updated by changing the function
    ``default_options`` attribute.

    .. [1] M. Nagahara and Y. Yamamoto, *Frequency-Domain Min-Max Optimization
       of Noise-Shaping Delta-Sigma Modulators*, IEEE Trans. SP, vol. 60 n. 6
       June 2012.

    For more information on the ``CVXOPT`` optimizer parameters, see
    the corresponding documentation. Similarly, for more infomration
    on the ``scs`` optimizer see its documentation. Finally, for more
    iniformation on the parameters controlling the optimization
    modeler and frontend ``cvxpy``, see its documentation.

    """
    # Manage optional parameters
    opts = digested_options(
        options, ntf_fir_minmax.default_options,
        ['show_progress', 'modeler'], [], False)
    dig_opts = {'show_progress': opts['show_progress'],
                'cvxpy_opts': {},
                'cvxpy_tdr_opts': {},
                'picos_opts': {}}
    if opts['modeler'] == 'cvxpy':
        opts.update(digested_options(
            options, ntf_fir_minmax.default_options,
            [], ['cvxpy_opts'], False))
        if opts['cvxpy_opts']['solver'] == 'cvxopt':
            dig_opts['cvxpy_opts'].update(digested_options(
                options, ntf_fir_minmax.default_options,
                [], ['cvxopt_opts'], False)['cvxopt_opts'])
            if opts['cvxpy_opts'].get('override_kktsolver', True):
                dig_opts['cvxpy_opts']['kktsolver'] = 'chol'
        elif opts['cvxpy_opts']['solver'] == 'scs':
            dig_opts['cvxpy_opts'].update(digested_options(
                options, ntf_fir_minmax.default_options,
                [], ['scs_opts'], False)['scs_opts'])
        opts['cvxpy_opts'].pop('override_kktsolver')
        dig_opts['cvxpy_opts'].update(opts['cvxpy_opts'])
        from ._fir_minmax_cvxpy import (
            ntf_fir_from_digested as _ntf_fir_from_digested)
    elif opts['modeler'] == 'cvxpy_old':
        dig_opts['cvxpy_tdr_opts'].update(digested_options(
            options, ntf_fir_minmax.default_options,
            [], ['cvxopt_opts'], False)['cvxopt_opts'])
        from ._fir_minmax_cvxpy_tdr import (
            ntf_fir_from_digested as _ntf_fir_from_digested)
    elif opts['modeler'] == 'picos':
        dig_opts['picos_opts'].update(digested_options(
            options, ntf_fir_minmax.default_options,
            [], ['cvxopt_opts'], False)['cvxopt_opts'])
        from ._fir_minmax_picos import (
            ntf_fir_from_digested as _ntf_fir_from_digested)
    else:
        raise ValueError('Unsupported modeling backend {}'.format(
            opts['modeler']))
    digested_options(options, {})
    if np.isscalar(f0):
        f0 = [f0]
    if np.isscalar(osr):
        osr = len(f0) * [osr]
    if not len(osr) == len(f0):
        raise ValueError('Incorrect multiband specification')
    # Do the computation
    ntf_ir = _ntf_fir_from_digested(order, osr, H_inf, f0, zf, **dig_opts)
    return (np.roots(ntf_ir), np.zeros(order), 1.)

ntf_fir_minmax.default_options = {"cvxpy_opts": {'override_kktsolver': False,
                                                 'solver': 'cvxopt'},
                                  "cvxopt_opts": {'maxiters': 100,
                                                  'abstol': 1e-7,
                                                  'reltol': 1e-6,
                                                  'feastol': 1e-6},
                                  'scs_opts': {'max_iters': 25000,
                                               'eps': 1e-10,
                                               'alpha': 1.8,
                                               'normalize': True,
                                               'use_indirect': False},
                                  'show_progress': True,
                                  'modeler': 'cvxpy_old'}


# Following part is deprecated


def synthesize_ntf_minmax(order=32, osr=32, H_inf=1.5, f0=0, zf=False,
                          **options):
    """
    Alias of :func:`ntf_fir_minmax`

    .. deprecated:: 0.11.0
       Function is now available from the :mod:`pydsm.NTFdesign`
       module with name :func:`ntf_fir_minmax`

    """
    warn("Function superseded by ntf_fir_minmax in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_fir_minmax(order, osr, H_inf, f0, zf, **options)
