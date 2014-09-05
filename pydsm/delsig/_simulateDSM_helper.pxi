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

# Helper inline functions for cython simulateDSM code

cdef inline double dbl_sat(double x, double a, double b):
    return a if x <= a else b if x>=b else x

cdef inline void ds_quantize(int N, double* y, int y_stride, \
    int* n, int n_stride, \
    double* v, int v_stride):
    """Quantize a signal according to a given number of levels."""
    cdef int qi
    cdef double L
    for qi in range(N):
        if n[qi*n_stride] % 2 == 0:
            v[qi*v_stride] = 2*floor(0.5*y[qi*y_stride])+1
        else:
            v[qi*v_stride] = 2*floor(0.5*y[qi*y_stride]+1)
        L = n[qi*n_stride]-1
        v[qi*v_stride]=dbl_sat(v[qi*v_stride],-L,L)

cdef inline void track_vabsmax(int N,\
    double* vabsmax, int vabsmax_stride,\
    double* x, int x_stride):
    cdef int i
    cdef double absx
    for i in range(N):
        absx=fabs(x[i*x_stride])
        if absx > vabsmax[i*vabsmax_stride]:
            vabsmax[i*vabsmax_stride]=absx

cdef inline double *dbldata(np.ndarray arr):
    return <double *>np.PyArray_DATA(arr)

cdef inline int *intdata(np.ndarray arr):
    return <int *>np.PyArray_DATA(arr)