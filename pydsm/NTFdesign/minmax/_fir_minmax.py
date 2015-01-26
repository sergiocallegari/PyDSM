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
    u"""
    Synthesize FIR NTF for LP or BP ΔΣ modulator by min-max optimization.

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

    Other parameters
    ----------------
    show_progress : bool, optional
        Provide extended output.
    modeler : string, optional
        modeling backend for the optimization problem. Currently, the
        ``cvxpy_old`` and the ``picos`` backends are supported.
        Default is ``cvxpy_old``.
    cvxopt_opts : dictionary, optional
        A dictionary of options for the ``cvxopt`` optimizer
        Allowed options include:

        ``maxiters``
            Maximum number of iterations (defaults to 100)
        ``abstol``
            Absolute accuracy (defaults to 1e-7)
        ``reltol``
            Relative accuracy (defaults to 1e-6)
        ``feastol``
            Tolerance for feasibility conditions (defaults to 1e-6)

        Do not use other options since they could break ``cvxpy`` in
        unexpected ways.

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

    See also
    --------
    cvxopt : for the optimizer parameters.
    """
    # Manage optional parameters
    opts = digested_options(
        options, ntf_fir_minmax.default_options,
        ['show_progress', 'modeler'], [], False)
    dig_opts = {'show_progress': opts['show_progress'],
                'cvxpy_opts': {},
                'tinoco_opts': {},
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
        dig_opts['tinoco_opts'].update(digested_options(
            options, ntf_fir_minmax.default_options,
            [], ['cvxopt_opts'], False)['cvxopt_opts'])
        from ._fir_minmax_tinoco import (
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
                                  'scs_opts': {'max_iters': 2500,
                                               'eps': 1e-12,
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
       Function is now available from the :mod:`NTFdesign` module with
       name :func:`ntf_fir_minmax`
    """
    warn("Function superseded by ntf_fir_minmax in "
         "NTFdesign module", PyDsmDeprecationWarning)
    return ntf_fir_minmax(order, osr, H_inf, f0, zf, **options)
