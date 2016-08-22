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

	def move(self):

		theta = 2.0*np.pi*np.random.rand(self.N)
		self.x += self.l*np.cos(theta)
		self.y += self.l*np.sin(theta)

class TimeIntegrator(object):
    
    def __init__(self,Brownian,Nt):
        self.N = Brownian.N
        self.Nt = Nt
        self.xData = np.zeros((Brownian.N,Nt),dtype=Brownian.dtype)
        self.yData = np.zeros((Brownian.N,Nt),dtype=Brownian.dtype)
        
    def forward(self,Brownian):
        for i in range(1,self.Nt):
            Brownian.move()
            self.xData[:,i]=Brownian.x
            self.yData[:,i]=Brownian.y
            
class QoI(object):
    
    def __init__(self):
        self.J=0.0
        
    def compute_J(self,TimeIntegrator):
        self.J = np.mean( np.square(TimeIntegrator.xData)+np.square(TimeIntegrator.yData) )