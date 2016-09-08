import numpy as np
import scipy as np

class Brownian(object):

    dtype = np.double

    def __init__(self,N,l,dt,
                        x0=None,y0=None):

        self.N = N
        self.l = l
        self.dt = dt

        self.x, self.y = x0, y0
        if x0 is None:
            self.x = np.zeros(N, dtype=self.dtype)
        if y0 is None:
            self.y = np.zeros(N, dtype=self.dtype)
            
        self.theta = np.zeros(N, dtype=self.dtype)
        
    def getTheta(self):
        self.theta = 2.0*np.pi*np.random.rand(self.N)        

    def move(self):
        #Note that l itself is the distance during dt, not the velocity
        self.x += self.l*np.cos(self.theta)
        self.y += self.l*np.sin(self.theta)

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
        
class TimeIntegrator(object):
    
    def __init__(self,Brownian,Nt):
        self.N = Brownian.N
        self.Nt = Nt
        self.dtype = Brownian.dtype
        
        self.xData = np.zeros((Nt,Brownian.N),dtype=self.dtype)
        self.yData = np.zeros((Nt,Brownian.N),dtype=self.dtype)
        
        self.theta = np.zeros((Nt,Brownian.N),dtype=self.dtype)
        
    def forward(self,Brownian):
        np.random.seed(3217)
        for i in range(self.Nt):
            Brownian.getTheta()
            self.theta[i,:] = Brownian.theta
            Brownian.move()
            self.xData[i,:]=Brownian.x
            self.yData[i,:]=Brownian.y
            
    def backward(self,AdjBrownian,QoI):
        
        if type(AdjBrownian)==particle_adj:
            QoI.dJdA = 0.0
            for i in range(self.Nt):
                nk = self.Nt-1-i
                QoI.compute_dJ_ptc(self,nk)
                AdjBrownian.x += AdjBrownian.dt*QoI.dx
                AdjBrownian.y += AdjBrownian.dt*QoI.dy
                QoI.compute_dJdA_ptc(AdjBrownian,self.theta[nk,:])
                
        elif type(AdjBrownian)==adj_particle:
            QoI.dJdA = 0.0
            for i in range(self.Nt):
                nk = self.Nt-1-i
                AdjBrownian.getTheta()
                AdjBrownian.move()
                QoI.compute_dJ_pdf(nk)
                AdjBrownian.N += len(QoI.dx)
                AdjBrownian.spwt = np.append(AdjBrownian.spwt,QoI.spwt)
                AdjBrownian.x = np.append(AdjBrownian.x,QoI.dx)
                AdjBrownian.y = np.append(AdjBrownian.y,QoI.dy)
                QoI.compute_dJdA_pdf(AdjBrownian,self,nk)
            
class QoI(object):
    
    def __init__(self,TimeIntegrator):
        self.dtype = TimeIntegrator.dtype
        self.N = TimeIntegrator.N
        self.Nt = TimeIntegrator.Nt
        self.J = 0.0
        self.dJdA = 0.0
#        self.dx = np.zeros(self.N,dtype=self.dtype)
#        self.dy = np.zeros(self.N,dtype=self.dtype)
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
        
class particle_adj(object):
    
    def __init__(self,Brownian):
        self.dtype = Brownian.dtype
        self.dt = -Brownian.dt
        self.N = Brownian.N
        self.x = np.zeros(self.N, dtype=self.dtype)
        self.y = np.zeros(self.N, dtype=self.dtype)
              
class adj_particle(Brownian):
    
    def __init__(self,ref_Br):
        Brownian.__init__(self,N=0,l=ref_Br.l,dt=ref_Br.dt,x0=[],y0=[])
        self.spwt = []