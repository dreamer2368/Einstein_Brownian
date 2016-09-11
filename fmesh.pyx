import numpy as np
cimport numpy as np

ctypedef np.float64_t DOUBLE

cdef extern:
    void c_ptc2pdf(int *Ng, double *dx, int *N, double *pg,
                    double *x, double *y, double *spwt_input)
    void c_pdfgrad(int *Ng, double *dx, double *pg,
                    double *pg_dx, double *pg_dy)

class mesh(object):
    
    dtype = np.double
    
    def __init__(self,Ng,L):
        self.Ng = Ng
        self.L = L
        self.dx = L/Ng
        self.xg = np.linspace(-L,L,2*Ng+1)
        self.yg = np.linspace(-L,L,2*Ng+1)

    def ptc2pdf(self,
                double[:] x, double[:] y,
                double[:] spwt=None):
        cdef int N, Ng
        cdef double dx
        N = len(x)
        Ng = self.Ng
        dx = self.dx
        cdef np.ndarray[DOUBLE,ndim=2] pg = np.zeros((2*self.Ng+1,2*self.Ng+1),dtype=np.double)
        if spwt is None:
            spwt = 1./N*np.ones(N)
        c_ptc2pdf(&Ng,&dx,&N,&pg[0,0],&x[0],&y[0],&spwt[0])
        return pg
    
    def pdfgrad(self,double[:,:] pg):
        cdef int Ng
        cdef double dx
        Ng = self.Ng
        dx = self.dx
        cdef np.ndarray[DOUBLE,ndim=2] pg_dx = np.zeros((2*Ng+1,2*Ng+1),dtype=np.double)
        cdef np.ndarray[DOUBLE,ndim=2] pg_dy = np.zeros((2*Ng+1,2*Ng+1),dtype=np.double)
        c_pdfgrad(&Ng,&dx,&pg[0,0],&pg_dx[0,0],&pg_dy[0,0])
        return pg_dx, pg_dy