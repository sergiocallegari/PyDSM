# -*- coding: utf-8 -*-

# Copyright (c) 2012, Sergio Callegari
# All rights reserved.
#
# Portions of code ported from the DELSIG toolbox
# Copyright (c) 2009, Richard Schreier

"""
DELSIG helpers routines
=======================
"""

import numpy as np
from math import ceil, sqrt
from _padding import padl, padr
from _tf import rmsGain
from ..utilities import db

__all__=["ds_synNTFobj1", "ds_f1f2", "ds_optzeros"]

def ds_synNTFobj1(x, p, osr, f0):
    """
    Objective function for synthesizeNTF.
    """
    z = np.exp(2j*np.pi*(f0+0.5/osr*x))
    if f0 > 0:
        z = padl(z,len(p)/2, np.exp(2j*np.pi*f0))
    z = np.concatenate((z, z.conj()))
    if f0 == 0:
        z = padr(z, len(p), 1)
    f1, f2 = ds_f1f2(osr, f0)

    return db(rmsGain((z, p, 1), f1, f2))

def ds_f1f2(OSR=64, f0=0, complex_flag=False):
    """
    Helper function.
    """
    if complex_flag:
        f1 = f0-0.5/OSR
        f2 = f0+0.5/OSR
    else:
        if f0 > 0.25/OSR:
            f1 = f0-0.25/OSR
            f2 = f0+0.25/OSR
        else:
            f1 = 0;
            f2 = 0.5/OSR
    return f1, f2
    
def ds_optzeros(n, opt=1):
    """
    Helper function for synthesizeNTF.
    
    Returns the zeros which minimize the in-band noise power of 
    a delta-sigma modulator's NTF.
    """
    if opt == 0:
        optZeros = np.zeros(ceil(n/2.))
    else:
        if n == 1:
            optZeros = 0.
        elif n == 2:
            if opt == 1:
                optZeros = sqrt(1./3)
            else:
                optZeros = 0.
        elif n == 3:
            optZeros = np.asarray([sqrt(3./5), 0.])
        elif n == 4:
            if opt == 1:
                discr = sqrt(9./49-3./35)
                tmp = 3./7
                optZeros = np.sqrt([tmp+discr, tmp-discr])
            else:
                optZeros = np.asarray([0., sqrt(5./7)])
        elif n == 5:
            discr = sqrt(25./81-5./21)
            tmp = 5./9
            optZeros = np.sqrt([tmp+discr, tmp-discr, 0.])
        elif n == 6:
            if opt == 1:
                optZeros = np.asarray([0.23862059, 0.66120988, 0.9324696])
            else:
                discr = sqrt(56.)/33
                tmp = 7./11
                optZeros = np.sqrt([0, tmp+discr, tmp-discr])
        elif n == 7:
            optZeros = np.asarray([0, 0.40584371, 0.74153078, 0.94910785])
        elif n == 8:
            if opt == 1:
                optZeros = np.asarray([0.18343709, 0.52553345, 0.79666684, 
                                       0.96028993])
            else:
                optZeros = np.asarray([0, 0.50563161, 0.79017286, 0.95914731])
        elif n == 9:
            optZeros = np.asarray([0, 0.32425101, 0.61337056, 0.83603082, 
                                   0.9681602])
        elif n == 10:
            if opt == 1:
                optZeros = np.asarray([0.1834370913, 0.5255334458, 
                                       0.7966668433, 0.9602899327])
            else:
                optZeros = np.asarray([0, 0.41572267, 0.67208682, 0.86238894,
                                       0.97342121])
        elif n == 11:
            optZeros = np.asarray([0, 0.26953955, 0.51909468, 0.73015137,
                                   0.88706238, 0.97822864])
        elif n == 12:
            if opt == 1:
                optZeros = np.asarray([0.12523875, 0.36783403, 0.58731921,
                                       0.7699033, 0.90411753, 0.9815607])
            else:
                optZeros = np.asarray([0, 0.35222363, 0.58006251, 0.76647993,
                                       0.90281326, 0.98132047])
        elif n == 13:
            optZeros = np.asarray([0, 0.23045331, 0.44849063, 0.64234828,
                                   0.8015776, 0.91759824, 0.98418306])
        elif n == 14:
            if opt == 1:
                optZeros = np.asarray([0.10806212, 0.31911586, 0.51525046,
                                       0.68729392, 0.82720185, 0.92843513,
                                       0.98628389])
            else:
                optZeros = np.asarray([0, 0.30524384, 0.50836649, 0.6836066,
                                       0.82537239, 0.92772336, 0.98615167])
        else:
            raise ValueError('Optimized zeros for n>14 are not available.')

    # Sort the zeros and replicate them.
    z = np.sort_complex(optZeros).reshape(-1)
    optZeros = np.zeros(n, complex)
    m = 0
    if (n%2) == 1:
        optZeros[0] = z[0]
        z = z[1:]
        m = m+1
    for i in xrange(len(z)):
        optZeros[m]   =  z[i]
        optZeros[m+1] = -z[i]
        m = m+2

    return optZeros
