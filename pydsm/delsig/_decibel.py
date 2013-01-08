# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.
#
# Portions of code ported from the DELSIG toolbox
# Copyright (c) 2009, Richard Schreier

"""
Collection of DELSIG style decibel routines
===========================================
"""

import numpy as np

__all__=["dbv", "dbp", "dbm", "undbv", "undbp", "undbm"]

def dbv(x):
    """
    Converts a voltage (or current) ratio to dB.

    Parameters
    ----------
    x : real or array_like of reals
        argument should be strictly positive

    Returns
    -------
    y : real or array_like of reals
        20*log10(x)
    """
    return 20*np.log10(x)

def dbp(x):
    """
    Converts a power ratio to dB.

    Parameters
    ----------
    x : real or array_like of reals
        argument should be strictly positive

    Returns
    -------
    y : real or array_like of reals
        10*log10(x)
    """
    return 10*np.log10(x)

def dbm(v, R=50.):
    """
    Converts argument from rms voltage to dBm.

    Parameters
    ----------
    v : real or array_like of reals
        rms value to be converted in dBm
    R : real, optional
        value of test resistor (defaults to 50)

    Returns
    -------
    y : real or array_like of reals
        dBm value corresponding to the power provided by the rms voltag v
        on the test resistor R, normalized over the reference value of 1mW
    """
    return 10*np.log10(np.abs(v**2)/R)+30

def undbv(x):
    """
    Converts argument from db to a voltage ratio.

    Inverse of `dbv`.

    Parameters
    ----------
    x : real or array_like of reals
        value to convert

    Returns
    -------
    y : real or array_like of reals
        10^(x/20)
    """
    return 10**(np.asarray(x)/20.)

def undbp(x):
    """
    Converts argument from db to a power ratio.

    Inverse of `dbp`.

    Parameters
    ----------
    x : real or array_like of reals
        value to convert

    Returns
    -------
    y : real or array_like of reals
        10^(x/10)
    """
    return 10**(np.asarray(x)/10.)

def undbm(x, R=50.):
    """
    Converts argument from dBm power to rms voltage.

    Inverse of `dbm`.

    Parameters
    ----------
    x : real or array_like of reals
        value to convert
    R : real, optional
        test resistor

    Returns
    -------
    y : real or array_like of reals
        rms voltage that applied over a resistor R gives a power equal to x
        in dBm (i.e. referred to 1mW)
    """
    return np.sqrt(R*10**((x-30)/10.))
