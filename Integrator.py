import numpy as np
import scipy as np
from Brownian import *
from mesh import *
from QoI import *
from adjoint import *
        
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