import numpy as np
import scipy as np
from Brownian import *
from adjoint import *
from Integrator import *
from fmesh import *

class QoI(object):
    
    def __init__(self,N,Nt):
        self.dtype = np.double
        self.N = N
        self.Nt = Nt
        self.J = 0.0
        self.dJdA = 0.0
        self.dx = None
        self.dy = None
        self.spwt = None
        
    def compute_J(self,TimeIntegrator):
        self.J = np.mean( np.square(TimeIntegrator.xData)+np.square(TimeIntegrator.yData) )
        
    def compute_dJ_ptc(self,TimeIntegrator,n):
        self.dx = 2.0*TimeIntegrator.xData[n,:]/self.Nt/self.N
        self.dy = 2.0*TimeIntegrator.yData[n,:]/self.Nt/self.N
        
    def compute_dJ_pdf(self,n):
        #radius of window function
        W = 3.5
        #mass increase during dt
        Z = -np.pi*W**4.0/2./self.Nt
        #uniform, constant specific weight
        Np = self.N*2/self.Nt
        self.spwt = Z/Np*np.ones(Np)
        r = np.random.rand(Np)
        r = np.power(r,1./4.)
        r = W*r
        theta = np.random.rand(Np)
        theta = 2.*np.pi*theta
        self.dx = r*np.cos(theta)
        self.dy = r*np.sin(theta)
        
    def compute_dJdA_ptc(self,AdjBrownian,theta):
        xFluct = np.cos(theta)
        yFluct = np.sin(theta)
        self.dJdA += ( np.inner(AdjBrownian.x,xFluct) + np.inner(AdjBrownian.y,yFluct) )/AdjBrownian.dt
        
    def compute_dJdA_pdf(self,AdjBrownian,TimeIntegrator,n):
        W = 6.
        Ng = 32
        dx = W/Ng
        qmesh = mesh(Ng,W)
        #original/adjoint pdf
        pg = qmesh.ptc2pdf(TimeIntegrator.xData[n,:],TimeIntegrator.yData[n,:])
        pg_hat = qmesh.ptc2pdf(AdjBrownian.x,AdjBrownian.y,AdjBrownian.spwt)
        #Gradients of pg
        pg_dx, pg_dy = qmesh.pdfgrad(pg)
        #Gradients of pg_hat
        pg_hat_dx, pg_hat_dy = qmesh.pdfgrad(pg_hat)
        self.dJdA += AdjBrownian.l/2.*np.sum( pg_dx*pg_hat_dx + pg_dy*pg_hat_dy )

def dJ_pdf_vw(self,n):
    #radius of window function
    W = 3.5
    #mass increase during dt
    Z = -np.pi*W**4.0/2./self.Nt
    #time-varying number of particles
    Np = np.floor( self.N/self.Nt*(1.+2.*n/self.Nt) )
    #uniform pdf
    r = np.random.rand(Np)
    r = np.power(r,1./2.)
    r = W*r
    theta = np.random.rand(Np)
    theta = 2.*np.pi*theta
    self.dx = r*np.cos(theta)
    self.dy = r*np.sin(theta)
    #time-varying, non-uniform specific weight
    self.spwt = Z/Np*2.*( np.square(self.dx)+np.square(self.dy) )/W/W