import numpy as np
import scipy as np

class Brownian(object):

	dtype = np.double

	def __init__(self,N,l,dt,
						x0=None,y0=None):

		self.N = N
		self.l = l
		self.dt = dt

		self.x0, self.y0 = x0, y0
		if x0 is None:
			self.x0 = np.zeros(N, dtype=self.dtype)
		if y0 is None:
			slef.y0 = np.zeros(N, dtype=self.dtype)

	def move(self):

		theta = 2.0*np.pi*np.random.rand(N)
		x0 += self.l*np.cos(theta)
		y0 += self.l*np.sin(theta)
