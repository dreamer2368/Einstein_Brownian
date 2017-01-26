import numpy as np
cimport numpy as np

ctypedef np.float64_t DOUBLE

cdef extern:
    void c_ptc2pdf(int *Ng, double *dx, double *dy, int *N, double *pg,
                    double *x, double *y, double *spwt_input)
    void c_ptc2pdf_tsc(int *Ng, double *dx, double *dy, int *N, double *pg,
                       double *x, double *y, double *spwt_input)
    void c_pdfgrad(int *Ng, double *dx, double *dy, double *pg,
                    double *pg_dx, double *pg_dy)
    void c_pdf2ptc(int *Ng, double *dx, double *dy, int *N, double *pg,
                   double *x, double *y, double *pp)
    void c_pdf2ptc_tsc(int *Ng, double *dx, double *dy, int *N, double *pg,
                       double *x, double *y, double *pp)

class mesh(object):

    dtype = np.double
    
    def __init__(self,Ng,Lx,Ly):
        self.Ng = Ng
        self.Lx = Lx
        self.Ly = Ly
        self.dx = Lx/Ng
        self.dy = Ly/Ng
        self.xg = np.linspace(-Lx,Lx,2*Ng+1)
        self.yg = np.linspace(-Ly,Ly,2*Ng+1)

    def ptc2pdf(self,
                double[:] x, double[:] y,
                double[:] spwt=None):
        cdef int N, Ng
        N = len(x)
        Ng = self.Ng
        cdef double dx,dy
        dx,dy = self.dx, self.dy
        if spwt is None:
            spwt = 1./N*np.ones(N)
        cdef double[:] sc=np.ascontiguousarray(spwt,dtype=np.double)
        cdef double[:] xc=np.ascontiguousarray(x,dtype=np.double)
        cdef double[:] yc=np.ascontiguousarray(y,dtype=np.double)
        pg = np.zeros((2*self.Ng+1,2*self.Ng+1))
        cdef double[:,:] pgc = np.ascontiguousarray(pg,dtype=np.double)
       
        c_ptc2pdf(&Ng,&dx,&dy,&N,&pgc[0,0],&xc[0],&yc[0],&sc[0])
        return pgc
    
    def ptc2pdf_tsc(self,
                    double[:] x, double[:] y,
                    double[:] spwt=None):
        cdef int N, Ng
        N = len(x)
        Ng = self.Ng
        cdef double dx,dy
        dx,dy = self.dx, self.dy
        if spwt is None:
            spwt = 1./N*np.ones(N)
        cdef double[:] sc=np.ascontiguousarray(spwt,dtype=np.double)
        cdef double[:] xc=np.ascontiguousarray(x,dtype=np.double)
        cdef double[:] yc=np.ascontiguousarray(y,dtype=np.double)
        pg = np.zeros((2*self.Ng+1,2*self.Ng+1))
        cdef double[:,:] pgc = np.ascontiguousarray(pg,dtype=np.double)
       
        c_ptc2pdf_tsc(&Ng,&dx,&dy,&N,&pgc[0,0],&xc[0],&yc[0],&sc[0])
        return pgc
    
    def pdf2ptc(self, double[:,:] pg,
                double[:] x, double[:] y):
        cdef int N, Ng
        N = len(x)
        Ng = self.Ng
        cdef double dx, dy
        dx, dy = self.dx, self.dy
        cdef double[:,:] pgc = np.ascontiguousarray(pg,dtype=np.double)
        cdef double[:] xc = np.ascontiguousarray(x,dtype=np.double)
        cdef double[:] yc = np.ascontiguousarray(y,dtype=np.double)
        cdef np.ndarray[DOUBLE,ndim=1] pp = np.zeros(N,dtype=np.double)
        cdef double[:] ppc = np.ascontiguousarray(pp,dtype=np.double)

        c_pdf2ptc(&Ng,&dx,&dy,&N,&pgc[0,0],&xc[0],&yc[0],&ppc[0])
        return ppc

    def pdf2ptc_tsc(self, double[:,:] pg,
                    double[:] x, double[:] y):
        cdef int N, Ng
        N = len(x)
        Ng = self.Ng
        cdef double dx, dy
        dx, dy = self.dx, self.dy
        cdef double[:,:] pgc = np.ascontiguousarray(pg,dtype=np.double)
        cdef double[:] xc = np.ascontiguousarray(x,dtype=np.double)
        cdef double[:] yc = np.ascontiguousarray(y,dtype=np.double)
        cdef np.ndarray[DOUBLE,ndim=1] pp = np.zeros(N,dtype=np.double)
        cdef double[:] ppc = np.ascontiguousarray(pp,dtype=np.double)

        c_pdf2ptc_tsc(&Ng,&dx,&dy,&N,&pgc[0,0],&xc[0],&yc[0],&ppc[0])
        return ppc

    def pdfgrad(self,double[:,:] pg):
        cdef int Ng
        Ng = self.Ng
        cdef double dx, dy
        dx, dy = self.dx, self.dy
        cdef double[:,:] pgc = np.ascontiguousarray(pg,dtype=np.double)
        pg_dx = np.zeros((2*Ng+1,2*Ng+1),dtype=np.double)
        pg_dy = np.zeros((2*Ng+1,2*Ng+1),dtype=np.double)
        cdef double[:,:] pg_dxc = np.ascontiguousarray(pg_dx,dtype=np.double)
        cdef double[:,:] pg_dyc = np.ascontiguousarray(pg_dy,dtype=np.double)

        c_pdfgrad(&Ng,&dx,&dy,&pgc[0,0],&pg_dxc[0,0],&pg_dyc[0,0])
        return pg_dxc, pg_dyc