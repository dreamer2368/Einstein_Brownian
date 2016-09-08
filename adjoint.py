import numpy as np
import scipy as np
from Brownian import *

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