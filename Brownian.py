import numpy as np
import scipy as np

class particle(object):
    dtype = np.double
    
    def __init__(self,N,dt,
                 x0=None,y0=None):
        self.N = N
        self.dt = dt
        self.x, self.y = x0, y0
        if x0 is None:
            self.x = np.zeros(N, dtype=self.dtype)
        if y0 is None:
            self.y = np.zeros(N, dtype=self.dtype)

class Brownian(particle):
    dtype = np.double

    def __init__(self,N,l,dt,
                        x0=None,y0=None):
        particle.__init__(self,N,dt,x0,y0)
        self.l = l
        self.theta = np.zeros(N, dtype=self.dtype)
        
    def getTheta(self):
        self.theta = 2.0*np.pi*np.random.rand(self.N)        

    def move(self):
        #Note that l itself is the distance during dt, not the velocity
        self.x += self.l*np.cos(self.theta)
        self.y += self.l*np.sin(self.theta)