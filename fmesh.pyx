import numpy as np
cimport numpy as np

ctypedef np.float64_t DOUBLE

cdef extern:
    void c_ptc2pdf(int *Ng, double *dx, int *N, double *pg,
                    double *x, double *y, double *spwt_input)
    void c_pdfgrad(int *Ng, double *dx, double *pg,
                    double *pg_dx, double *pg_dy)

def f_ptc2pdf(int Ng,
                double dx,
                int N,
                double[:] x, double[:] y, double[:] spwt_input):
    cdef np.ndarray[DOUBLE,ndim=2] pg = np.zeros((2*Ng+1,2*Ng+1),dtype=np.double)
    c_ptc2pdf(&Ng,&dx,&N,&pg[0,0],&x[0],&y[0],&spwt_input[0])
    return pg

def f_pdfgrad(int Ng, double dx,
                double[:,:] pg):
    cdef np.ndarray[DOUBLE,ndim=2] pg_dx = np.zeros((2*Ng+1,2*Ng+1),dtype=np.double)
    cdef np.ndarray[DOUBLE,ndim=2] pg_dy = np.zeros((2*Ng+1,2*Ng+1),dtype=np.double)
    c_pdfgrad(&Ng,&dx,&pg[0,0],&pg_dx[0,0],&pg_dy[0,0])
    return pg_dx, pg_dy

class mesh(object):
    
    dtype = np.double
    
    def __init__(self,Ng,L):
        self.Ng = Ng
        self.L = L
        self.dx = L/Ng
        self.xg = np.linspace(-L,L,2*Ng+1)
        self.yg = np.linspace(-L,L,2*Ng+1)

    def ptc2pdf(self,x,y,spwt=None):
        N = len(x)
        if spwt is None:
            spwt = 1./N*np.ones(N)
        pg = f_ptc2pdf(self.Ng,self.dx,N,x,y,spwt)
        return pg
    
    def pdfgrad(self,pg):
        pg_dx, pg_dy = f_pdfgrad(self.Ng,self.dx,pg)
        return pg_dx, pg_dy