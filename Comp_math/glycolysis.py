# Number X.9.10

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
 
class glyc:
	solution = []
	def __init__(self, alpha, beta, y_0, t_0, T_k):
		self.alpha = alpha
		self.beta = beta
		self.y_0 = y_0
		self.segm = [t_0, T_k]

	def F(self, t, y):
		F = [1 - y[0] * y[1]]
		F.append(self.alpha * y[1] * (y[0] - (1 + self.beta) / (y[1] + self.beta)))
		return np.array(F)

	def scheme(self, tau):
		equations = lambda y: np.array(y) - self.solution[-1]- tau * (self.F(0, y) + self.F(0, self.solution[-1])) / 2
		return equations

	def non_lyn(self, tau):
		return fsolve(self.scheme(tau), self.solution[-1])

	def get_solution(self, N):
		tau = (self.segm[1] - self.segm[0]) / N
		self.solution.append(np.array(self.y_0))
		for i in range(1, N):
			self.solution.append(np.array(self.non_lyn(tau)))
		return self.solution

T_k = 50
l = glyc(1000, 10, [1, 1e-3], 0, T_k)
N = 10000
q = l.get_solution(N)
plt.plot(1, 1e-3, 'ro')
plt.plot([i[0] for i in q], [i[1] for i in q])
plt.show()
