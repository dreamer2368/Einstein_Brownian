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
        QoI.dJdA = 0.0
        for i in range(self.Nt):
            nk = self.Nt-1-i
            QoI.compute_dJ(self,nk)
            AdjBrownian.x += AdjBrownian.dt*QoI.dx
            AdjBrownian.y += AdjBrownian.dt*QoI.dy
            QoI.compute_dJdA(AdjBrownian,self.theta[nk,:])
            
class QoI(object):
    
    def __init__(self,TimeIntegrator):
        self.dtype = TimeIntegrator.dtype
        self.N = TimeIntegrator.N
        self.Nt = TimeIntegrator.Nt
        self.J = 0.0
        self.dJdA = 0.0
        self.dx = np.zeros(self.N,dtype=self.dtype)
        self.dy = np.zeros(self.N,dtype=self.dtype)
        
    def compute_J(self,TimeIntegrator):
        self.J = np.mean( np.square(TimeIntegrator.xData)+np.square(TimeIntegrator.yData) )
        
    def compute_dJ(self,TimeIntegrator,n):
        self.dx = 2.0*TimeIntegrator.xData[n,:]/self.Nt/self.N
        self.dy = 2.0*TimeIntegrator.yData[n,:]/self.Nt/self.N
        
    def compute_dJdA(self,AdjBrownian,theta):
        xFluct = np.cos(theta)
        yFluct = np.sin(theta)
        self.dJdA += ( np.inner(AdjBrownian.x,xFluct) + np.inner(AdjBrownian.y,yFluct) )/AdjBrownian.dt
        
class AdjBrownian(object):

    def __init__(self,Brownian):
        
        self.dtype = Brownian.dtype
        self.N = Brownian.N
        self.dt = -Brownian.dt

        self.x = np.zeros(self.N, dtype=self.dtype)
        self.y = np.zeros(self.N, dtype=self.dtype)