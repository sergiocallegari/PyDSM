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
