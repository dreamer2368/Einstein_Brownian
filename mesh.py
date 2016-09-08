import numpy as np
import scipy as np

class mesh(object):
    
    dtype = np.double
    
    def __init__(self,Ng,L):
        self.Ng = Ng
        self.L = L
        self.dx = L/Ng
        self.xg = np.linspace(-L,L,2*Ng+1)
        self.yg = np.linspace(-L,L,2*Ng+1)

    def ptc2pdf(self,x,y,spwt=None):
        pg = np.zeros((2*self.Ng+1,2*self.Ng+1),dtype=self.dtype)
        N = len(x)
        if spwt is None:
            spwt = 1./N*np.ones(N)
        for i in range(N):
            xp = x[i]
            yp = y[i]
            xgl = np.floor(xp/self.dx)
            xgr = xgl + 1
            ygl = np.floor(yp/self.dx)
            ygr = ygl + 1
            pg[xgl+self.Ng,ygl+self.Ng] += spwt[i]*(1.-xp+xgl*self.dx)*(1.-yp+ygl*self.dx)
            pg[xgl+self.Ng,ygr+self.Ng] += spwt[i]*(1.-xp+xgl*self.dx)*(yp-ygl*self.dx)
            pg[xgr+self.Ng,ygl+self.Ng] += spwt[i]*(xp-xgl*self.dx)*(1.-yp+ygl*self.dx)
            pg[xgr+self.Ng,ygr+self.Ng] += spwt[i]*(xp-xgl*self.dx)*(yp-ygl*self.dx)
        pg *= 1./self.dx
        return pg
    
    def pdfgrad(self,pg):
        pg_dx = np.zeros((2*self.Ng+1,2*self.Ng+1),dtype=self.dtype)
        pg_dy = np.zeros((2*self.Ng+1,2*self.Ng+1),dtype=self.dtype)
        pg_dx[1:2*self.Ng,:] = ( pg[2:2*self.Ng+1,:]-pg[0:2*self.Ng-1,:] )/2./self.dx
        pg_dy[:,1:2*self.Ng] = ( pg[:,2:2*self.Ng+1]-pg[:,0:2*self.Ng-1] )/2./self.dx
        return pg_dx, pg_dy